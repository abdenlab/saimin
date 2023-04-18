use pyo3::prelude::*;
use pyo3::types::PyBytes;

use std::str;
use std::sync::Arc;

use arrow::array::{ArrayRef, Int32Array, StringArray};
use arrow::ipc::writer::FileWriter;
use arrow::record_batch::RecordBatch;

#[pyclass]
struct BamReader {
    reader: Option<bam::IndexedReader<std::fs::File>>,
}

#[pymethods]
impl BamReader {
    #[new]
    fn py_new(path: &str) -> Self {
        let reader = bam::IndexedReader::from_path(path).unwrap();
        Self { reader: Some(reader) }
    }

    pub fn fetch(&mut self, chrom: &str, start: u32, end: u32) -> PyObject {
        let reader = self.reader.as_mut().expect("file closed.");
        let ref_id = reader
            .header()
            .reference_id(chrom)
            .expect("Invalid reference name");
        let region = bam::Region::new(ref_id, start, end);
        // let names = reader.header().reference_names();
        let region_viewer = reader.fetch(&region).unwrap();
        Python::with_gil(|py| {
            let py_bytes = PyBytes::new(py, &write_ipc(region_viewer)).into();
            py_bytes
        })
    }

    /// https://github.com/PyO3/pyo3/issues/1205#issuecomment-1164096251 for advice on `__enter__`
    pub fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }

    pub fn close(&mut self) {
        self.reader = None;
    }

    pub fn __exit__(&mut self, _exc_type: PyObject, _exc_value: PyObject, _traceback: PyObject) {
        self.close();
    }
}

fn write_ipc(
    region_viewer: bam::bam_reader::RegionViewer<std::fs::File>,
    // reference_names: &[String],
) -> Vec<u8> {
    let mut ipc = Vec::new();

    // TODO: get ref_id as category
    let mut ref_id: Vec<i32> = vec![];
    let mut start: Vec<i32> = vec![];
    let mut end: Vec<i32> = vec![];
    let mut name: Vec<Option<String>> = vec![];
    let mut cigar: Vec<Option<String>> = vec![];
    let mut seq: Vec<Option<String>> = vec![];
    let mut qual: Vec<Option<String>> = vec![];

    // Map row-wise entries to column-wise arrays
    for record in region_viewer {
        let record = record.unwrap();
        ref_id.push(record.ref_id());
        start.push(record.start());
        end.push(record.calculate_end());
        name.push(Some(str::from_utf8(record.name()).unwrap().to_string()));
        cigar.push(Some(record.cigar().to_string()));
        seq.push(Some(str::from_utf8(&record.sequence().to_vec()).unwrap().to_string()));
        qual.push(Some(str::from_utf8(&record.qualities().raw()).unwrap().to_string()));
    }

    // https://docs.rs/bam/latest/bam/record/struct.Record.html
    let batch = RecordBatch::try_from_iter(vec![
        ("ref_id", Arc::new(Int32Array::from(ref_id)) as ArrayRef),
        ("start", Arc::new(Int32Array::from(start)) as ArrayRef),
        ("end", Arc::new(Int32Array::from(end)) as ArrayRef),
        ("name", Arc::new(StringArray::from(name)) as ArrayRef),
        ("cigar", Arc::new(StringArray::from(cigar)) as ArrayRef),
        ("seq", Arc::new(StringArray::from(seq)) as ArrayRef),
        ("qual", Arc::new(StringArray::from(qual)) as ArrayRef),
    ])
    .unwrap();

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
