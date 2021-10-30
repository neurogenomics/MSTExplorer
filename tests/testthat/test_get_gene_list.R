geneData <- data.frame(
  "list_name" = c("a","a","a","a","b","b","b","c","c","c","b","b","c","c","c","d","d"),
  "Gene" = c("x","y","z","x","x","x","y","y","z","x","x","y","y","z","x","x","y"),
  "other_col" = seq(1,17)
)

test_that("returns correct genes",{
  expect_equal(get_gene_list("a",geneData, "list_name","Gene"),geneData[geneData$list_name == "a","Gene"])
  expect_equal(get_gene_list("b",geneData, "list_name","Gene"),geneData[geneData$list_name == "b","Gene"])
  expect_equal(get_gene_list("c",geneData, "list_name","Gene"),geneData[geneData$list_name == "c","Gene"])
  expect_equal(get_gene_list("d",geneData, "list_name","Gene"),geneData[geneData$list_name == "d","Gene"])
  expect_equal(get_gene_list("e",geneData, "list_name","Gene"),geneData[geneData$list_name == "e","Gene"], options(warn = -1))
  expect_warning(get_gene_list("e",geneData, "list_name","Gene"),paste("gene list", "e", "is not present in" ,"list_name"))

})

test_that("error when nonexistant columns are selected", {
  expect_error(get_gene_list("d",geneData, "Wrongcol","Gene"), "undefined columns selected")
  expect_error(get_gene_list("d",geneData, "list_name","Wrongcol"), "undefined columns selected")
  expect_error(get_gene_list("d",geneData, "Wrongcol","Wrongcol"), "undefined columns selected")
})
