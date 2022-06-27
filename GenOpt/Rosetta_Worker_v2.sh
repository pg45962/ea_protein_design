#!/bin/bash -e

### Updated to use the most recent ROSETTA build, as well as take glycosilations
### into account. -- LPARA 16102020

target_dir="$1"                            
if [ -z "$target_dir" ]                 
then
    echo -e "No target directory provided."
    exit 1
fi
target_seq="$2"                        
if [ -z "$target_seq" ]                
then
    echo -e "No sequence provided."
    exit 1
fi

FILE="design.xml ab_rdb.pdb"
output_name="$target_seq"
n_replicas=5
PATH_TO_COMP=/compilations/rosetta_bin_linux_2020.08.61146_bundle/main/source/bin
EXE=rosetta_scripts.mpi.linuxgccrelease
minims_dir="$target_dir/minims"
final_out_dir="$minims_dir/$target_seq"

### Initialize working directory.
work_dir=$(mktemp -d -t ci-XXXXXXXXXX --tmpdir="$minims_dir")
cp -p "$target_dir"/const_param/design.xml "$work_dir/"         
cp -p "$target_dir"/const_param/ab_rdb.pdb "$work_dir/"       
./Edit_Design.py -i "$work_dir"/design.xml -s "$target_seq"   

### Run rosetta
cd "$work_dir"
mpirun -np 1 -quiet "$PATH_TO_COMP"/"$EXE" -s "$FILE" -parser:protocol design.xml -nstruct "$n_replicas" -out:file:silent "$output_name"_structures.silent -out:file:scorefile "$output_name"_score.sc -use_input_sc -ex1 -ex2 -ex2aro -database /compilations/rosetta_bin_linux_2020.08.61146_bundle/main/database -include_sugars -alternate_3_letter_codes pdb_sugar -maintain_links -auto_detect_glycan_connections -load_PDB_components false> file.log 2>&1 || true 
energy=$(sort -nk 7 "$output_name"_score.sc | head -n 1 | awk '{print "Score:", $3, "  I_sc:", $7}')
isc=$(sort -nk 7 "$output_name"_score.sc | head -n 1 | awk '{print $7}')
best_pdb=$(sort -nk 7 "$output_name"_score.sc | head -n 1 | awk '{print $NF}')
"$PATH_TO_COMP"/extract_pdbs.mpi.linuxgccrelease -in::file::silent "$output_name"_structures.silent -in::file::tags "$best_pdb" -include_sugars -alternate_3_letter_codes pdb_sugar -maintain_links -auto_detect_glycan_connections > file_extract.log 2>&1
cd - > /dev/null

if [ ! -f "$work_dir"/"$output_name"_score.sc ]
then
    #rm -rf "$work_dir"
    exit 0
fi

### Check if final output directory exists
if [ ! -f "$final_out_dir"/lowest_score.txt ]        ### Se a pasta com os mins nao existir
then
    rm -rf "$final_out_dir"
    mkdir "$final_out_dir"
    echo "$energy" > "$final_out_dir"/lowest_score.txt
    cp "$work_dir"/"$best_pdb".pdb "$final_out_dir"/ab_rdb.pdb
    cat "$final_out_dir"/lowest_score.txt
    rm -rf "$work_dir"
    exit 0
fi

### If it exists compare the new isc to the previous best
oldisc=$(cat "$final_out_dir"/lowest_score.txt | awk '{print $4}')
if (( $(echo "$isc < $oldisc" | bc -l) ))   
then
    echo "$energy" > "$final_out_dir"/lowest_score.txt
    cp "$work_dir"/"$best_pdb".pdb "$final_out_dir"/ab_rdb.pdb
    cat "$final_out_dir"/lowest_score.txt
    rm -rf "$work_dir"
    exit 0
else
    cat "$final_out_dir"/lowest_score.txt
    rm -rf "$work_dir"
    exit 0
fi
