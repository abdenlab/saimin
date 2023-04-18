use pyo3::prelude::*;
use pyo3::types::PyBytes;

use std::str;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, GenericBinaryBuilder, GenericStringArray, Int32Array, StringDictionaryBuilder,
};
use arrow::datatypes::Int8Type;
use arrow::ipc::writer::FileWriter;
use arrow::record_batch::RecordBatch;

use noodles::bam;
use noodles::bgzf;
use noodles::sam;

#[pyclass]
struct BamReader {
    reader: bam::IndexedReader<bgzf::Reader<std::fs::File>>,
    header: sam::Header,
}

#[pymethods]
impl BamReader {
    #[new]
    fn py_new(path: &str) -> Self {
        let mut reader = bam::indexed_reader::Builder::default()
            .build_from_path(path)
            .unwrap();
        let header = reader.read_header().unwrap();
        Self { reader, header }
    }

    pub fn fetch(&mut self, chrom: &str, start: u32, end: u32) -> PyObject {
        let region = format!("{}:{}-{}", chrom, start, end).parse().unwrap();
        let query = self.reader.query(&self.header, &region).unwrap();
        let ipc = write_ipc(query);
        Python::with_gil(|py| PyBytes::new(py, &ipc).into())
    }

    /// https://github.com/PyO3/pyo3/issues/1205#issuecomment-1164096251 for advice on `__enter__`
    pub fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }

    pub fn __exit__(&mut self, _exc_type: PyObject, _exc_value: PyObject, _traceback: PyObject) {}
}

fn write_ipc(query: bam::reader::Query<std::fs::File>) -> Vec<u8> {
    let mut refs = StringDictionaryBuilder::<Int8Type>::new();

    let size = 100; // TODO: get size from region_viewer?
    let mut starts = Int32Array::builder(size);
    let mut ends = Int32Array::builder(size);

    let mut names = GenericBinaryBuilder::<i32>::new();
    let mut cigars = GenericBinaryBuilder::<i32>::new();
    let mut seqs = GenericBinaryBuilder::<i32>::new();
    let mut quals = GenericBinaryBuilder::<i32>::new();

    // Map row-wise entries to column-wise arrays
    for result in query {
        let record = result.unwrap();
        // TODO: map reference id to name
        refs.append("chr1").unwrap();
        starts.append_value(record.alignment_start().unwrap().get() as i32);
        ends.append_value(record.alignment_end().unwrap().get() as i32);
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
fn bram(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<BamReader>()?;
    Ok(())
}
