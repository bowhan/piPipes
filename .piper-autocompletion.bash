# TAB auto-completion for pipipe
_piper () {
	local cur
	COMPREPLY=( "install" "small" "small2" "rnaseq" "rnaseq2" "deg" "chip" "chip2" "genome" ) 
  return 0
}
complete -F _piper  ./piper

