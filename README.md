Scripts to download and process data from the pridedb dtabases (ionbot results).

How to use:

run the queries from "beekeeper-queries" in a psql shell to download the data.
Gzip the results and then run the python scripts in order:

* 01: map peptides to obtain position of PTMs on the proteins; generate a file with the mapped peptidoforms and their respective psm counts.
* 02: calculate absolute PTM counts
* 03: calculate relative PTM counts (#PSMs with modified residue / total coverage of the residue)
* 04: remove duplicated psms from PSM file
* 05: Add flashLFQ quantification to PSMs
* 06: Combine files from 05 into one (make a single scrip??)
* 07: 
* 08:
* 09: