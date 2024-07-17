use std::cmp::max;

use jseqio::writer::SeqRecordWriter;
use rand::{distributions::Distribution, RngCore, SeedableRng};

fn sample_reads(db: &jseqio::seq_db::SeqDB, read_length: usize, n_reads: usize, seed: u64) -> Vec<Vec<u8>> {

    // Create a seeded RNG
    let mut seed_32_bytes: [u8; 32] = [0; 32];
    seed_32_bytes[0..8].copy_from_slice(&seed.to_le_bytes());
    let mut rng = rand::rngs::StdRng::from_seed(seed_32_bytes);

    // We want to pick a sequence with probability proportional to the number of (read_length)-mers in the sequence.
    let seq_weights: Vec<usize> = db.iter().map(|rec| 
        max(rec.seq.len() as isize - read_length as isize + 1, 0) as usize)
    .collect();
    let seq_distribution = rand::distributions::WeightedIndex::new(seq_weights).unwrap();

    let mut sampled_reads = Vec::<Vec::<u8>>::with_capacity(n_reads);
    while sampled_reads.len() < n_reads {
        let seq_idx = seq_distribution.sample(&mut rng);
        let seq = db.get(seq_idx).seq.to_vec();

        assert!(seq.len() >= read_length); // This should have weight zero and thus never happen

        let start = rng.next_u64() as usize % (seq.len() - read_length + 1);
        let read = seq[start..start+read_length].to_owned();

        if read.iter().any(|&c| c != b'A' && c != b'C' && c != b'G' && c != b'T') {
            continue;
        }

        sampled_reads.push(read);
    }

    sampled_reads

}



fn main() {

    let cli = clap::Command::new("Sample k-mers. k-mers with non-ACGT characters are not sampled.")
        .arg_required_else_help(true)
        .arg(clap::Arg::new("source-file")
            .long("source-file")
            .short('s')
            .value_parser(clap::value_parser!(std::path::PathBuf))
        )
        .arg(clap::Arg::new("output-file")
            .long("output-file")
            .short('o')
            .value_parser(clap::value_parser!(std::path::PathBuf))
        )
        .arg(clap::Arg::new("k")
            .short('k')
            .value_parser(clap::value_parser!(usize))
        )
        .arg(clap::Arg::new("howmany")
            .long("howmany")
            .short('n')
            .value_parser(clap::value_parser!(usize))
        );

    let matches = cli.get_matches();
    let source_file = matches.get_one::<std::path::PathBuf>("source-file").unwrap();
    let output_file = matches.get_one::<std::path::PathBuf>("output-file").unwrap();
    let k = *matches.get_one::<usize>("k").unwrap();
    let howmany = *matches.get_one::<usize>("howmany").unwrap();

    let reader = jseqio::reader::DynamicFastXReader::from_file(source_file).unwrap();
    let seq_db = reader.into_db().unwrap();


    // Random seed from current time
    let seed = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_micros();
    let seed = (seed % u64::MAX as u128) as u64;

    let reads = sample_reads(&seq_db, k, howmany, seed);

    let mut writer = jseqio::writer::DynamicFastXWriter::new_to_file(output_file).unwrap();

    for seq in reads {
        writer.write_ref_record(&jseqio::record::RefRecord{head: b"", seq: &seq, qual: None}).unwrap();
    }

}
