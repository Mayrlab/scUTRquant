## Download Rules
##
## Any rules added here will automatically be imported into the main Snakefile.
## Important notes:
##   - all `input`/`output` must be relative to main Snakefile
##   - `conda` is relative to this file
##   - `shell` will execute relative to main Snakefile
##   - `download_10X_whitelists` shows how to run script in-place (see shell)

rule download_10X_whitelists:
    input: "extdata/bxs/download_10X_whitelists.sh"
    output:
        "extdata/bxs/737K-april-2014_rc.txt",
        "extdata/bxs/737K-august-2016.txt",
        "extdata/bxs/3M-february-2018.txt"
    conda: "../../envs/downloading.yaml"
    shell:
        """
        pushd $(dirname {input})
        sh $(basename {input})
        popd
        """
