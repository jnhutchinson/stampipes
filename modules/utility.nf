/// This file is only for "utility" processes that are extremely generic.

process publish {

  publishDir params.outdir

  input:
    tuple val(filename), path(infile)
  output:
    path(filename)

  script:
    """
    ln -s "$infile" "$filename"
    """
}
