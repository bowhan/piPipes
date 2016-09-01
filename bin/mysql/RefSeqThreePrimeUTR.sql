SELECT
  chrom,
  IF(strand = '+', cdsEnd, txStart) AS start,
  IF(strand = '+', txEnd, cdsStart) AS end,
  name,
  score,
  strand
FROM refGene;