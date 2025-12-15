# Data Directory

This directory is for storing input data files.

## Important Notes

- **Do not commit large data files** to the repository
- Large data files (>.fits, .csv, .dat, .hdf5) are excluded via .gitignore
- Consider using Git LFS for large files or hosting them elsewhere

## Data Structure

Suggested organization:
```
data/
├── raw/              # Original, unprocessed data
├── processed/        # Cleaned and preprocessed data
└── sample/           # Small sample datasets for testing
```

## Obtaining Data

Add instructions here for:
- Where to download the data
- How to preprocess it
- Expected format and structure

## Data Format

Document the expected format of input data:
- File types (FITS, CSV, etc.)
- Required columns/fields
- Units and conventions
