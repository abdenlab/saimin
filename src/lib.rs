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

#[pyclass]
struct BamReader {
    reader: bam::IndexedReader<std::fs::File>,
}

#[pymethods]
impl BamReader {
    #[new]
    fn py_new(path: &str) -> Self {
        let reader = bam::IndexedReader::from_path(path).unwrap();
        Self { reader }
    }

    pub fn fetch(&mut self, chrom: &str, start: u32, end: u32) -> PyObject {
        let header = self.reader.header().clone();
        let ref_id = header.reference_id(chrom).expect("Invalid reference name.");
        let region = bam::Region::new(ref_id, start, end);
        let ipc = write_ipc(self.reader.fetch(&region).unwrap(), &header);
        Python::with_gil(|py| PyBytes::new(py, &ipc).into())
    }

    /// https://github.com/PyO3/pyo3/issues/1205#issuecomment-1164096251 for advice on `__enter__`
    pub fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }

    pub fn __exit__(&mut self, _exc_type: PyObject, _exc_value: PyObject, _traceback: PyObject) {}
}

fn write_ipc(
    region_viewer: bam::bam_reader::RegionViewer<std::fs::File>,
    header: &bam::Header,
) -> Vec<u8> {
    let mut refs = StringDictionaryBuilder::<Int8Type>::new();

    let size = 100; // TODO: get size from region_viewer?
    let mut starts = Int32Array::builder(size);
    let mut ends = Int32Array::builder(size);

    let mut names = GenericBinaryBuilder::<i32>::new();
    let mut cigars = GenericBinaryBuilder::<i32>::new();
    let mut seqs = GenericBinaryBuilder::<i32>::new();
    let mut quals = GenericBinaryBuilder::<i32>::new();

    // Map row-wise entries to column-wise arrays
    for record in region_viewer {
        let record = record.unwrap();
        refs.append(header.reference_name(record.ref_id() as u32).unwrap())
            .unwrap();
        starts.append_value(record.start());
        ends.append_value(record.calculate_end());
        names.append_value(record.name());
        cigars.append_value(record.cigar().to_string());
        seqs.append_value(record.sequence().to_vec());
        quals.append_value(record.qualities().raw());
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
