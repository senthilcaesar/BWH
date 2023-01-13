export CLICOLOR=1
export LSCOLORS=GxFxCxDxBxegedabagaced
alias eris='ssh -Y -X sq566@eristwo.partners.org'
FFTW=/opt/homebrew/Cellar/fftw/3.3.10_1
LGBM_PATH=/Users/sq566/Downloads/tmp/LightGBM
#alias luna='/Users/sq566/Downloads/luna-base-0.27/luna'




# Check another version /usr/local/bin/luna
# Github token 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
# docker run -d -p 3838:3838 nsrr-cohort-matrix
# rsync -r -v --progress -e ssh sq566@eristwo.partners.org:/data/nsrr/working/apples-edfs/\*.edf . ( remote to local )
# glpat-S78KuPeS7acxSuNErYVx ( gitlab-nsrr-project-token )
# scp sq566@eristwo.partners.org:/data/nsrr/working/apples-edfs/\*.edf
# sshpass -p 'SonyXperia$9' rsync -r -v --progress -e ssh sq566@eristwo.partners.org:/data/nsrr/datasets/ccshs/polysomnography/edfs/ccshs-trec-1800786.edf .

alias moon='ssh -i "/Users/sq566/Downloads/test_rshinyserver.pem" ubuntu@ec2-18-188-74-28.us-east-2.compute.amazonaws.com' #( NAP Server )
# ssh -i "remnrem_website.pem" ubuntu@ec2-3-23-150-12.us-east-2.compute.amazonaws.com ( NAP Web app )



# bash -c "nohup sh gilla.sh &" ( Run jobs in background )
# Copy local files to EC2 instance
# scp -i test_rshinyserver.pem ml.zip ubuntu@ec2-18-188-74-28.us-east-2.compute.amazonaws.com:
# cmake -DCMAKE_CXX_COMPILER=g++-12 -DCMAKE_C_COMPILER=gcc-12 .. ( For LGBM )
#
#
#
#
#
# Update moonlight code to EC2
# scp -i ~/Downloads/test_rshinyserver.pem -r moonlight ubuntu@ec2-18-188-74-28.us-east-2.compute.amazonaws.com:~/.
