# Command to build conda package
```bash
conda mambabuild conda.recipe/
conda convert --platform osx-64 </path/to/package.tar.bz2> -o <outputdir/>
anaconda upload <converted_package>
```