params.scriptfile = ""

if (params.scriptfile == "") {
  die
}

scripts = Channel.create()
file(params.scriptfile).readLines().each { scripts << it }
scripts.close()

process run_script {
  input:
  val script_name from scripts

  errorStrategy "ignore"

  script:
  """
  cd "\$(dirname "${script_name}")"
  bash "${script_name}"
  """
}
