## Convert Rd to md
## Link: https://rdrr.io/cran/Rd2md/man/Rd2markdown.html

## set directory to folder with files
setwd("/Users/adrianasistig/Documents/GitHub/Giotto/rd_md")

filelist = list.files(pattern = ".Rd")
  for (i in 1:length(filelist)) {
    input <- filelist[i]
    #print(input)
    output <- str_replace(input, ".Rd",".md")
    #output = sub(pattern = "\\.Rd$", replacement = ".md", x = i)
    #print(output)
    Rd2markdown(input, output)
  }

## to convert md to rst
# in unix
## go into directory with files and run:
# m2r2 *.md
# rm *.md

## Everything should be an .rst now
# You'll need the m2r2 extenion on your local to run it (and it must be added to conf.py)
