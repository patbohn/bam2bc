use std::fs::{File};
use std::io::{BufReader, BufRead, Write};
use std::path::PathBuf;
use std::str::from_utf8;
//use std::io::{BufRead, BufReader};
//use flate2::read::MultiGzDecoder;
use std::error::Error;
use std::convert::TryInto;

use flate2::write::GzEncoder;
use flate2::Compression;

use structopt::StructOpt;

use bam::RecordReader;
use bam::record::Record;
use bam::record::tags::{TagValue};

use stringmetrics::levenshtein_limit;

///get phred (intergers) to error probability table
pub mod phred_int_to_prob;
use phred_int_to_prob::PHRED_TO_ERROR_PROB;

#[derive(Debug, StructOpt)]
#[structopt(name = "bampA", about = "Detects and estimates polyA/T tail lengths in DNA sequencing data with a given adapter sequence.")]
struct Config {
    #[structopt(parse(from_os_str), short="i", help = "Input bam file")]
    input_file: PathBuf,

    #[structopt(short="o", long, help = "Output statistics file", default_value ="output")]
    output_file: String,
    
    #[structopt(short="a", long, help = "file containing barcode sequences, one per row, now empty lines, no annotations")]
    barcode_file: String,
    
    #[structopt(short="c", long, help = "constant part of sequence that is located 5' of the variable barcode", default_value="CATGGTCATAGCTGTTTCCTG")]
    constant_part: String,
    
    #[structopt(short="m", long, help = "number of tolerated mismatches to find the constant sequence", default_value="10")]
    num_mismatches: u8,
    
    #[structopt(short="n", long, help = "number of nucleotides to search the adapter within at the ends", default_value="200")]
    adapter_search_region: usize,
    
    #[structopt(short="e", long, help = "number of nucleotides to ignore at the 3' end to calculate mean qscore of read (ignore accidental DNA part basecalls)", default_value="90")]
    ignore_end_for_qscore: usize,
    
    #[structopt(short="s", long, help = "Ignore moves with samples above this threshold to calculate mean translocation speed (correcting for stalling/tails)", default_value="200")]
    max_samples_per_move: usize,
    
    #[structopt(short="v", long, help = "Will print to stdout")]
    verbose: bool,
}

const THREADS: u16 = 8;

fn main() -> std::io::Result<()> {
    
    let config = Config::from_args();
    
    //define output paths
    let per_read_filename = format!("{}.tsv.gz", config.output_file);
    
    //create output files
    let per_read_statsfile = File::create(per_read_filename)?;
    let mut writer = GzEncoder::new(per_read_statsfile, Compression::fast());
    
    //read in bam file
    let mut bam_reader = bam::BamReader::from_path(config.input_file, THREADS).unwrap();
    //read bam let reader = BufReader::new(MultiGzDecoder::new(input_file));
    
    // read in config settings
    let constant_part: String = config.constant_part;
    let max_acceptable_dist: u8 = config.num_mismatches;
    let adapter_search_region: usize = config.adapter_search_region;
    let ignore_end_for_qscore: usize = config.ignore_end_for_qscore;
    let max_samples_per_move: usize = config.max_samples_per_move;
    let verbose: bool = config.verbose;

    // read in barcodes to search for
    let barcodes = read_in_barcodes(&config.barcode_file).expect("Failed to read barcode file");
    let num_barcodes = barcodes.len();
    
    //read_id, read_len, unique_adapter_detected, adapter_location, tail_type, tail_stride_length, nt_stride_speed, est_tail_length
    let barcode_mindist_cols: Vec<String> = (1..=num_barcodes).map(|i| format!("\tmindist_barcode{}", i)).collect();
    let barcode_loc_cols: Vec<String> = (1..=num_barcodes).map(|i| format!("\tloc_barcode{}", i)).collect();
    
    writeln!(&mut writer, "read_id\tread_len\tsamples_per_nt\ttransloc_stdev\tmean_phred\tmean_phred_RNA\tconst_dist\tconst_loc{}{}", barcode_mindist_cols.join(""), barcode_loc_cols.join(""))?;
    
    //initialize record and other variables here to save on allocation
    let mut record = bam::Record::new();
    
    loop {
        // reader: impl RecordReader
        // New record is saved into record.
        match bam_reader.read_into(&mut record) {
            // No more records to read.
            Ok(true) => {},
            Ok(false) => break,
            Err(e) => panic!("{}", e),
        }
        
        // Do somethind with the record.
        let read_id = from_utf8(record.name()).unwrap_or("NA");
        let raw_seq = record.sequence();
        let str_vec = raw_seq.to_vec();
        let read_seq = String::from_utf8_lossy(&str_vec);
        let read_len = read_seq.len();
        
        //let mut nt_transloc_speed: f32 = 0.0;
        //let mut nt_transloc_stdev: f32 = 0.0;
        let mut barcode_search_start = if read_len > adapter_search_region {read_len - adapter_search_region} else {0};
        
        if verbose {println!("read_id: {}\tread_seq: {}\tread_len: {}\tadapter_search_region: {}\tbarcode_search_start: {}", read_id, read_seq, read_len, adapter_search_region, barcode_search_start)}
        //if read is in mRNA orientation then adapter should be located at 3' end
        
        let (const_dist, const_loc) = search_for_constant(&read_seq, &constant_part, max_acceptable_dist);
        if const_dist < max_acceptable_dist {
            if verbose {println!("M13 barcode most likely detected at {} with distance {}. Adjusting start of search for barcode.", const_loc, const_dist)}
            barcode_search_start = const_loc;
            if read_len < (barcode_search_start + constant_part.len()){
                if verbose {println!("Read likely lacks 3' end")}
            }
        }
        
        let (bc_dists, bc_loc) = compare_w_barcodes(&read_seq[barcode_search_start..], &barcodes);
        
        let bc_dist_str = bc_dists.iter().map(|x| x.to_string()).collect::<Vec<String>>().join("\t");
        let bc_loc_str = bc_loc.iter().map(|x| (x+barcode_search_start).to_string()).collect::<Vec<String>>().join("\t");
        
        if verbose {println!("{}", bc_dist_str)}
        if verbose {println!("{}", bc_loc_str)}
        
        let moves_durations = get_nt_durations(&record, &verbose);
        
        let (nt_transloc_speed, nt_transloc_stdev) = calc_mean_move_per_nt(&moves_durations, max_samples_per_move);
        
        let quality = record.qualities().raw();
        
        let mean_error_prob = calc_mean_error(&quality);
        let mean_quality = error_prob_to_phred(mean_error_prob);
        
        let mut mean_quality_RNA = -1.0;
        
        if read_len > ignore_end_for_qscore {
            let DNA_start: usize = read_len - ignore_end_for_qscore;
            let mean_error_prob = calc_mean_error(&quality[..DNA_start]);
            mean_quality_RNA = error_prob_to_phred(mean_error_prob);
        } 
        
        
        writeln!(&mut writer, "{}\t{}\t{:.2}\t{:.2}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}", read_id, read_len, nt_transloc_speed, nt_transloc_stdev, mean_quality, mean_quality_RNA, const_dist, const_loc, bc_dist_str, bc_loc_str)?;
        }
    
    Ok(())
}


fn read_in_barcodes(barcode_file: &str) -> Result<Vec<String>, Box<dyn Error>> {
    
    let file = File::open(&barcode_file)?;
    let reader = BufReader::new(file);
    
    let mut barcodes = Vec::new();
    
    for line in reader.lines() {
        let line_contents = line?;
        barcodes.push(line_contents.to_uppercase().replace("U", "T"));
    }
    println!("Read in {} barcodes, assuming all are off length {}", barcodes.len(), barcodes[0].len());
    Ok(barcodes)
}

//compares sequence in window with all barcodes and returns lowest match
fn compare_w_barcodes(read_seq : &str, barcodes : &Vec<String>) -> (Vec<usize>, Vec<usize>) {
    let bc_len = barcodes[0].len();
    let mut bc_dists = vec![bc_len as usize; barcodes.len()];
    let mut bc_loc = vec![0 as usize; barcodes.len()];
    
    //extend seq to ensure we can also compare for 3' trimmed seq
    let mut extended_seq = read_seq.to_string();
    for _ in 0..bc_len {
        extended_seq.push('X');
    }
    
    let read_seq = &extended_seq;
    let last_slice: usize = read_seq.len() - bc_len;
    
    for slice_index in 0..last_slice {
        let read_slice = &read_seq[slice_index..(slice_index+bc_len)];
        for bc in 0..barcodes.len() {
            
            let lev_dist: usize = levenshtein_limit(&read_slice, &barcodes[bc], bc_len as u32).try_into().unwrap();
            if lev_dist < bc_dists[bc] {
                bc_dists[bc] = lev_dist;
                bc_loc[bc] = slice_index;
            }
        }
    }
    (bc_dists, bc_loc)
}

fn search_for_constant(read_seq : &str, constant_part : &String, max_acceptable_dist : u8) -> (u8, usize) {
    let bc_len = constant_part.len();
    let mut bc_dist: u8 = max_acceptable_dist;
    let mut bc_loc:usize = 0;
    
    //extend seq to ensure we can also compare for 3' trimmed seq
    let mut extended_seq = read_seq.to_string();
    for _ in 0..bc_len {
        extended_seq.push('X');
    }
    
    let read_seq = &extended_seq;
    let last_slice: usize = read_seq.len() - bc_len;
    
    for slice_index in 0..last_slice {
        let read_slice = &read_seq[slice_index..(slice_index+bc_len)];
        let lev_dist: u8 = levenshtein_limit(&read_slice, &constant_part, max_acceptable_dist as u32).try_into().unwrap();
        if lev_dist < bc_dist {
            bc_dist = lev_dist;
            bc_loc = slice_index;
        }
    }
    (bc_dist, bc_loc)
}


//https://github.com/nanoporetech/dorado/issues/337
//move table is found at tag "mv", it is int starting with "c,stride,move0..."
//we should be able to trim the first non-boolean value and then read in the rest as bool
//to estimate the duration of each base we can iterate through the move table, resetting count at index when a new 1 occurs

fn get_nt_durations(record: &Record, verbose: &bool) -> Vec<usize> {
    if let Some(tag_value) = record.tags().get(b"mv") {
        match tag_value {
            
            TagValue::IntArray(tag_array) => {
                if *verbose {println!("Detected mv table")}
                calc_moves_per_nt(tag_array.raw(), verbose)},
            _ => {panic!("Unexpected tag type at mv");}
        }
    } else {
        panic!("Could not read in tag!")
    }
}

fn calc_moves_per_nt(moves_table: &[u8], verbose: &bool) -> Vec<usize> {
    if *verbose {println!("calculating moves per nt")}
    let stride: usize = moves_table[0].try_into().unwrap();

    let mut stay_vector: Vec<usize> = Vec::new();
    let mut count: usize = 1;

    for mv in &moves_table[1..] {
        match mv {
            0 => {
                count +=1;
            }
            1 => {
                stay_vector.push(count*stride);
                count = 1;
            }
            _ => {
                panic!("Unexpected value in move table! Expected only values of <0,1>")
            }
        }
    }
    if *verbose {println!("Found {} moves", stay_vector.len())}
    stay_vector
}


//perform counting of each occurence to get frequencies
fn calc_mean_move_per_nt(duration_per_nt: &Vec<usize>, max_samples_per_move : usize) -> (f32,f32) {
    // Find the maximum and minimum values in the array
    let max_val = *duration_per_nt.iter().max().unwrap_or(&0);
    let min_val = 0;

    // Create a counting array to store the count of each element
    let mut count = vec![0; (max_val - min_val + 1) as usize];
    
    let mut sum_wo_outlier: f32 = 0.0;
    let mut count_wo_outlier: f32 = 0.0;
    // Count the occurrences of each element in the input array
    for &num in duration_per_nt.iter() {
        let index = (num) as usize;
        count[index] += 1;
        if num < max_samples_per_move {
            sum_wo_outlier += num as f32;
            count_wo_outlier += 1.0;
        }
    }
    let mean_move_wo_outlier: f32 = sum_wo_outlier / count_wo_outlier;
    
    let actual_max = if max_val > max_samples_per_move {max_samples_per_move} else {max_val};
    let sum_of_squared_diff: f32 = count[..actual_max as usize].iter().enumerate().map(|(i, &freq)| (i as f32 -mean_move_wo_outlier).powi(2) * freq as f32).sum();
    let stdev_wo_outlier: f32 = (sum_of_squared_diff / count_wo_outlier).sqrt();

    (mean_move_wo_outlier, stdev_wo_outlier)
}


fn calc_mean_error(quality_slice: &[u8]) -> f64 {
    
    let sum_error: f64 = quality_slice.iter().map(|&x| PHRED_TO_ERROR_PROB[x as usize]).sum();
    return sum_error / quality_slice.len() as f64
}

fn error_prob_to_phred(prob: f64) -> f64 {
    return -10.0_f64 * prob.log10()
}