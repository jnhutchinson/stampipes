/// This file is only for "utility" processes that are extremely generic.

process publish {

  publishDir params.outdir

  input:
    tuple val(filename), path("__infile__")
  output:
    path(filename)

  script:
    """
    ln -s "__infile__" "$filename"
    """
}
