# Full-Year Clonality Runbook

This runbook is for running a full year of clonality data, month by month, and then combining the month workbooks into one Excel file with an overview.

You can use the same flow for `2025` or `2024`.
Only update the paths and run names to match the machine you are on.

## 1. Run A Full Year

### Full 2025
```bash
python3 /path/to/OUS/scripts/run_clonality_2025_monthly.py \
  --input-root /path/to/data/Klonalitet/2025_data \
  --output-root /path/to/output/full_year_runs \
  --run-name full_2025_run \
  --max-workers 1 \
  --folder-workers 1
```

### Full 2024
```bash
python3 /path/to/OUS/scripts/run_clonality_2025_monthly.py \
  --input-root /path/to/data/Klonalitet/2024_data \
  --output-root /path/to/output/full_year_runs \
  --run-name full_2024_run \
  --max-workers 1 \
  --folder-workers 1
```

Notes:
- `max-workers 1` and `folder-workers 1` are the safest settings on laptops and offline machines with limited memory.
- The script still works for 2024 even though the filename says `2025`.
- Output is written into a fresh run folder, so older runs are preserved.

## 2. Resume If The Run Stops

### Resume 2025
```bash
python3 /path/to/OUS/scripts/run_clonality_2025_monthly.py \
  --input-root /path/to/data/Klonalitet/2025_data \
  --output-root /path/to/output/full_year_runs \
  --run-name full_2025_run \
  --resume-existing \
  --max-workers 1 \
  --folder-workers 1
```

### Resume 2024
```bash
python3 /path/to/OUS/scripts/run_clonality_2025_monthly.py \
  --input-root /path/to/data/Klonalitet/2024_data \
  --output-root /path/to/output/full_year_runs \
  --run-name full_2024_run \
  --resume-existing \
  --max-workers 1 \
  --folder-workers 1
```

## 3. Combine The Month Workbooks Into One Excel File

### Combine 2025
```bash
python3 /path/to/OUS/scripts/combine_clonality_2025_overview.py \
  --run-root /path/to/output/full_year_runs/full_2025_run \
  --year-label 2025
```

### Combine 2024
```bash
python3 /path/to/OUS/scripts/combine_clonality_2025_overview.py \
  --run-root /path/to/output/full_year_runs/full_2024_run \
  --year-label 2024
```

Default output:
- `track-clonality-2025-overview.xlsx`
- `track-clonality-2024-overview.xlsx`

You can also choose a custom output path:

```bash
python3 /path/to/OUS/scripts/combine_clonality_2025_overview.py \
  --run-root /path/to/output/full_year_runs/full_2025_run \
  --year-label 2025 \
  --output /path/to/output/full_year_runs/full_2025_run/my-full-2025-overview.xlsx
```

## 4. What The Combined Workbook Contains

The combined workbook includes:
- `Overview`
- `Patient_Runs_<YEAR>`
- `Control_Runs_<YEAR>`
- `PK_Peaks_<YEAR>`
- `Weak_Ladders_<YEAR>`
- `PK_Outliers_<YEAR>`

This gives:
- all patient rows for the year
- all control rows for the year
- all PK rows for the year
- a quick list of the remaining weak ladder cases
- a quick list of PK marker outliers

## 5. Important Output Folders

Inside a run folder such as `full_2025_run`, you will get:
- `full_2025_run_manifest.json`
- `month_folder_lists/`
- `month_runs/2025_01 ... 2025_12`

Inside each month folder:
- `track-clonality.xlsx`
- `backfill_state.json`
- `run_summary.json`
- `analysis/`
- `feature_artifacts/`
- `candidate_artifacts/`

## 6. Linux Offline Build

To build the offline Linux app bundle:

```bash
cd /path/to/OUS
./packaging/build_linux.sh
```

Expected artifact:
- `dist/releases/Fraggler_Linux_offline.zip`

Run it on the offline Linux PC:

```bash
unzip Fraggler_Linux_offline.zip -d ~/Fraggler
cd ~/Fraggler/Fraggler_Linux
chmod +x Fraggler
./Fraggler
```

## 7. Recommended Workflow On The Offline Linux PC

1. Build `Fraggler_Linux_offline.zip` on a machine with Docker and internet access.
2. Copy the zip to the offline Linux PC.
3. Extract and launch the app.
4. Copy the data folder to the Linux PC.
5. Run the full-year command from a terminal on that machine.
6. Combine the month workbooks into one overview workbook when the run is finished.
