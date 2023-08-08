# Command to build conda package
```bash
mamba mambabuild conda.recipe/
mamba convert --platform osx-64 </path/to/package.tar.bz2> -o <outputdir/>
anaconda upload <converted_package>
```