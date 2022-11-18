cd /Users/sq566

git clone https://github.com/remnrem/luna-base.git

cd luna-base

make ARCH=MAC FFTW=/opt/homebrew/Cellar/fftw/3.3.10_1 -j4 LGBM=1 LGBM_PATH=/Users/sq566/Downloads/tmp/LightGBM

echo alias luna='/Users/sq566/luna-base/luna' >> /Users/sq566/.bashrc
echo alias destrat='/Users/sq566/luna-base/destrat' >> /Users/sq566/.bashrc
echo alias behead='/Users/sq566/luna-base/behead' >> /Users/sq566/.bashrc

source /Users/sq566/.bashrc
