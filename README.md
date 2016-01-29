# rna-seq
* /build directory is ignored

* SeqTK added as a submodule.  Use these additional commands during a clone:
```
git submodule init
git submodule update
```

Optionally, you may clone recursively:
```
git clone --recursive https://github.com/UNR-CFB/rna-seq.git
```

INSTALLING ALL THE TOOLS

```
samtools
cd src/samtools
./configure
```
fix anything that screams at you
sudo yum ?????

```
make
make prefix=../../build/samtools/ install
cd ../../build/samtools/bin
echo 'export PATH=$PATH:'`pwd` >> $HOME/.bash_profile
> above is for jr-style path insertion
```

# Committing changes
The following commands assume you are operating from the project root directory.

Pull the newest changes from the repository before attempting to push:
```
git pull origin
```

Stage your next commit:
```
git add .
```

Check the staged commit contains all the expected changes:
```
git status
```

Commit the changes:
```
git commit
```

A text editor will start.  Input a description of your changes and save.

Push your changes to Github:
```
git push origin
```
