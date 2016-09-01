SELECT
  chrom,
  IF(strand = '+', txStart, cdsEnd) AS start,
  IF(strand = '+', cdsStart, txEnd) AS end,
  name,
  score,
  strand
FROM refGene;