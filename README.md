A command-line tool for designing, evaluating, and visualizing PCR primers in both singleplex and multiplex assays. It combines Primer3, thermodynamic QC, BLAST specificity checks (with local fallback), and rich output—ASCII maps, color-coded metrics, and optional graphical plots and Excel exports.
To Clone Repository for Python use...
    git clone https://github.com/joshuajenk/primerdiscovery.git
    cd primerdiscovery
    Install Repositories --> pip install pandas primer3-py biopython matplotlib thefuzz openpyxl
    Run Tool --> python main.py

To Run Tool...
s --> singleplex reaction, five optimal pairs for one qPCR/RT-PCR
m --> multiplex reaction, five optimal multiplex panels for up to 10+ pathogens.
templates --> templates can be added to the 'data' folder -- currently filled with 30+ blood-borne pathogens
custom --> add custom nucleic acid sequence, 100 to 1000 bp
ncbi --> ncbi accession codes can be found on NCBI PubMed per each gene... NOTE: entire gene will be searched, not specific targets

