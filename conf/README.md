## How to add new profile configurations

1. Choose a short (a few words) unique name for the profile that only has lowercase letters and underscores. This will be referred to as the `<profile ID>` here.
2. Add a new CSV with the input metadata to `test/data/metadata/`, named `<Profile ID>.csv`
3. Copy an existing profile file in `conf/`, rename it to `<profile ID>.config`, and replace any mentions of the copied profile with the profile ID within the text of the file.
4. Add a line to `nextflow.config` inside the `profiles {` section with the syntax `<Profile ID>   { includeConfig 'conf/<Profile ID>.config' }`
5. (Optional) If you think this profile would be useful to make test data sets for testing the main report outside of the pipeline, you can add the `<Profile ID>` to `scripts/make_report_test_datasets.sh`, in the line that looks like `PROFILES=(xanthomonas <Profile ID>)`.
