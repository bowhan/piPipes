# TAB auto-completion for pipipe
_pipipe () {
	local cur
	COMPREPLY=( "install small" "small2" "rnaseq" "rnaseq2" "deg" "chip" "chip2" "genome" ) 
  return 0
}
complete -F _pipipe  ./pipipe

