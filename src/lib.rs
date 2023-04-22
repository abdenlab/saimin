use pyo3::prelude::*;
use pyo3::types::PyBytes;

use std::str;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, GenericBinaryBuilder, GenericStringArray, Int32Array, StringDictionaryBuilder, UInt8Array,
};
use arrow::datatypes::Int8Type;
use arrow::ipc::writer::FileWriter;
use arrow::record_batch::RecordBatch;

use noodles::bam;
use noodles::bgzf;
use noodles::core::{Position, Region};
use noodles::sam;

type BufferedReader = std::io::BufReader<std::fs::File>;

#[pyclass]
struct BamReader {
    reader: bam::IndexedReader<bgzf::Reader<BufferedReader>>,
    header: sam::Header,
}

#[pymethods]
impl BamReader {
    #[new]
    fn py_new(path: &str) -> Self {
        let index = bam::bai::read(format!("{}.bai", path)).unwrap();
        let file = std::fs::File::open(path).unwrap();
        let bufreader = std::io::BufReader::with_capacity(1024 * 1024, file);
        let mut reader = bam::indexed_reader::Builder::default()
            .set_index(index)
            .build_from_reader(bufreader)
            .unwrap();
        let header = reader.read_header().unwrap();
        Self { reader, header }
    }

    pub fn fetch(&mut self, chrom: &str, start: usize, end: usize) -> PyObject {
        let start = Position::try_from(start + 1).unwrap();
        let end = Position::try_from(end).unwrap();
        let query = self
            .reader
            .query(&self.header, &Region::new(chrom, start..=end))
            .unwrap();
        let ipc = write_ipc(query, &self.header);
        Python::with_gil(|py| PyBytes::new(py, &ipc).into())
    }

    pub fn fetch_noodles(&mut self, chrom: &str, start: usize, end: usize) -> () {
        let start = Position::try_from(start + 1).unwrap();
        let end = Position::try_from(end).unwrap();
        let query = self
            .reader
            .query(&self.header, &Region::new(chrom, start..=end))
            .unwrap();
        for record in query {
            let _record = record.unwrap();
        }
        print!("done");
    }

    /// https://github.com/PyO3/pyo3/issues/1205#issuecomment-1164096251 for advice on `__enter__`
    pub fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }

    pub fn __exit__(&mut self, _exc_type: PyObject, _exc_value: PyObject, _traceback: PyObject) {}
}

fn write_ipc(query: bam::reader::Query<BufferedReader>, header: &sam::Header) -> Vec<u8> {
    let mut refs = StringDictionaryBuilder::<Int8Type>::new();

    let size = 100; // TODO: get size from region_viewer?
    let mut starts = Int32Array::builder(size);
    let mut ends = Int32Array::builder(size);
    let mut mapqs = UInt8Array::builder(size);
    let mut names = GenericBinaryBuilder::<i32>::new();
    let mut cigars = GenericBinaryBuilder::<i32>::new();
    let mut seqs = GenericBinaryBuilder::<i32>::new();
    let mut quals = GenericBinaryBuilder::<i32>::new();

    // Map row-wise entries to column-wise arrays
    for result in query {
        let record = result.unwrap();
        let ref_name = match record.reference_sequence(&header) {
            Some(Ok((name, _))) => name,
            None => "unknown",
            _ => panic!("error getting reference sequence"),
        };
        refs.append(ref_name).unwrap();
        starts.append_value(record.alignment_start().unwrap().get() as i32);
        ends.append_value(record.alignment_end().unwrap().get() as i32);
        mapqs.append_value(record.mapping_quality().unwrap().get());
        names.append_value(record.read_name().unwrap());
        cigars.append_value(record.cigar().to_string());
        seqs.append_value(record.sequence().to_string());
        quals.append_value(record.quality_scores().to_string());
    }

    let names = GenericStringArray::try_from_binary(names.finish()).unwrap();
    let cigars = GenericStringArray::try_from_binary(cigars.finish()).unwrap();
    let seqs = GenericStringArray::try_from_binary(seqs.finish()).unwrap();
    let quals = GenericStringArray::try_from_binary(quals.finish()).unwrap();

    // https://docs.rs/bam/latest/bam/record/struct.Record.html
    let batch = RecordBatch::try_from_iter(vec![
        ("ref", Arc::new(refs.finish()) as ArrayRef),
        ("start", Arc::new(starts.finish()) as ArrayRef),
        ("end", Arc::new(ends.finish()) as ArrayRef),
        ("mapq", Arc::new(mapqs.finish()) as ArrayRef),
        ("name", Arc::new(names) as ArrayRef),
        ("cigar", Arc::new(cigars) as ArrayRef),
        ("seq", Arc::new(seqs) as ArrayRef),
        ("qual", Arc::new(quals) as ArrayRef),
    ])
    .unwrap();

    let mut ipc = Vec::new();
    {
        let cursor = std::io::Cursor::new(&mut ipc);
        let mut writer = FileWriter::try_new(cursor, &batch.schema()).unwrap();
        writer.write(&batch).unwrap();
        writer.finish().unwrap();
    }
    ipc
}

/// A Python module implemented in Rust.
#[pymodule]
fn saimin(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<BamReader>()?;
    Ok(())
}
