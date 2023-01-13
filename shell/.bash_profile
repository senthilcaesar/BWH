
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/sq566/opt/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/sq566/opt/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/Users/sq566/opt/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/Users/sq566/opt/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

export PATH="/Users/sq566/.rbenv/versions/3.0.4/bin:$PATH"
export CLICOLOR=1
export LSCOLORS=GxFxCxDxBxegedabagaced


# Enable new docker build backend
export DOCKER_BUILDKIT=1
alias eris='ssh -Y -X sq566@eristwo.partners.org'
alias matlab='/Applications/MATLAB_R2022b.app/bin/matlab'
HISTFILESIZE=10000000
HISTSIZE=10000000
export PATH="/Applications/Sublime Text.app/Contents/SharedSupport/bin:$PATH"

alias eris='ssh -Y -X sq566@eristwo.partners.org'

alias luna='/Users/sq566/Programme/luna-base/luna'
alias destrat='/Users/sq566/Programme/luna-base/destrat'
alias behead='/Users/sq566/Programme/luna-base/behead'
alias fixrows='/Users/sq566/Programme/luna-base/fixrows'
