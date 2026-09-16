from pathlib import Path
import sys

import numpy as np
import pgenlib


root = Path(__file__).resolve().parents[3]
asset = root / "snp_spectral_foundation" / "pilot_assets" / "geuvadis_tensorqtl_chr18"
pgen = asset / "GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18.pgen"

with pgenlib.PgenReader(str(pgen).encode()) as reader:
    print({"samples": reader.get_raw_sample_ct(), "variants": reader.get_variant_ct()})
    values = np.empty((reader.get_raw_sample_ct(), 10), dtype=np.int8)
    reader.read_range(0, 10, values, sample_maj=1)
    print({"shape": values.shape, "values": sorted(np.unique(values).tolist())})

