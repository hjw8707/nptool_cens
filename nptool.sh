#!/bin/sh

# test if export is supported
export 1>/dev/null 2>/dev/null
if [ $? -eq 0 ]; then
  CMD="export"
  SEP="="
else
  setenv 1>/dev/null 2>/dev/null
  if [ "${?} == 0" ]; then
  CMD="setenv"
  SEP=" "
  else
  echo "Neither setenv nor export found!"
  fi
fi 

# find script path
if [ -n "$ZSH_VERSION" ]; then
   SCRIPTPATH="$( cd "$( dirname "${(%):-%x}" )" && pwd )"
elif [ -n "$tcsh" ]; then
   SCRIPTPATH="$( cd "$( dirname "$0" )" && pwd )"
elif [ -n "$BASH_VERSION" ]; then
   SCRIPTPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
else
   echo "neither bash or zsh is used, abort"
   exit 1
fi

# export NPTOOL environment variable
${CMD} NPTOOL${SEP}$SCRIPTPATH

NPARCH=$(uname)
# mac os x case
if [ "${NPARCH}" = "Darwin" ] ; 
then
  ${CMD} DYLD_LIBRARY_PATH${SEP}$NPTOOL/NPLib/lib:$DYLD_LIBRARY_PATH
  ${CMD} DYLD_LIBRARY_PATH${SEP}$NPTOOL/NPSimulation/lib:$DYLD_LIBRARY_PATH
else 
  ${CMD} LD_LIBRARY_PATH${SEP}$NPTOOL/NPLib/lib:$LD_LIBRARY_PATH
  ${CMD} LD_LIBRARY_PATH${SEP}$NPTOOL/NPSimulation/lib:$LD_LIBRARY_PATH
fi

${CMD} PATH=$NPTOOL/NPLib/bin:$PATH
${CMD} PATH=$NPTOOL/NPSimulation/bin:$PATH

alias npt='cd $NPTOOL'  
alias npl='cd $NPTOOL/NPLib'  
alias nps='cd $NPTOOL/NPSimulation'
${CMD} npa_not_supported='npa is now longer supported, use npp instead'
alias npa='echo $npa_not_supported'


# open a project
function npp {
  if [[ $1 == *"Example"* ]]
  then
    cd $NPTOOL/Examples/$1 
  else
    cd $NPTOOL/Projects/$1
  fi
}
# tab completion for npp
_npp() {
  # Pointer to current completion word.
  local cur
  # Array variable storing the possible completions.
  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"

  # Combine subdirectories from both Projects and Examples
  local proj_dirs example_dirs
  proj_dirs=$(basename -a "$NPTOOL"/Projects/* 2>/dev/null)
  example_dirs=$(basename -a "$NPTOOL"/Examples/* 2>/dev/null)

  local LIST="${proj_dirs} ${example_dirs}"

  COMPREPLY=( $(compgen -W "$LIST" -- "$cur") )

  return 0
}

# associate the tab completion to npp
if [ -n "$ZSH_VERSION" ]; then
  # Reuse the bash completion function in zsh.
  autoload -U +X bashcompinit && bashcompinit
  complete -F _npp -o filenames npp
else
  # the rest of the world use standard posix complete
  complete -F _npp -o filenames npp
fi

${CMD} Geant4_DIR${SEP}$G4LIB
${CMD} NPLib_DIR${SEP}$NPTOOL/NPLib
