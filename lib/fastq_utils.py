#!/usr/bin/env python

import copy
import gzip
import json
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.request as request
from contextlib import closing
from multiprocessing import Process
from pathlib import Path
import requests

from fqutil_api import authenticateByEnv, getHostManifest

# Default bowtie2 threads.
BT2_THREADS = 2
SAM_THREADS = 1

# hisat2 has problems with spaces in filenames
# prevent spaces in filenames. if one exists link the file to a no-space version.


def link_space(file_path):
    result = file_path
    name = os.path.splitext(os.path.basename(file_path))[0]
    if " " in name:
        clean_name = name.replace(" ", "")
        result = file_path.replace(name, clean_name)
        if not os.path.exists(result):
            os.symlink(file_path, result)
            # subprocess.run(["ln", "-s", file_path, result], check=True)
    return result


def get_param(program, feature, tool_params):
    if program in tool_params and feature in tool_params[program]:
        return tool_params[program][feature]
    return None


def run_fastqc(read_list, output_dir, job_data, tool_params):
    rcount = 0
    threads = get_param("fastqc", "-p", tool_params)
    fastqc_base_cmd = ["fastqc", "--memory", "4096", "-t", str(threads) if threads else "8"]
    for r in read_list:
        rcount += 1
        if len(r.get("fastqc", [])) == 0:
            fastqc_cmd = fastqc_base_cmd + ["--outdir", output_dir]
            if "read2" in r:
                fastqc_cmd += [r["read1"], r["read2"]]
                r["fastqc"].append(
                    os.path.join(
                        output_dir, os.path.basename(r["read1"]) + ".fastqc.html"
                    )
                )
                r["fastqc"].append(
                    os.path.join(
                        output_dir, os.path.basename(r["read2"]) + ".fastqc.html"
                    )
                )
            else:
                fastqc_cmd += [r["read1"]]
                r["fastqc"].append(
                    os.path.join(
                        output_dir, os.path.basename(r["read1"]) + ".fastqc.html"
                    )
                )
            print((" ".join(fastqc_cmd)))
            subprocess.run(fastqc_cmd, check=True)


def find_prefix(filename):
    r1_parts = filename.split(".")
    prefix_pos = (
        -2
        if (
            (r1_parts[-1] == "gz" or r1_parts[-1] == "gzip")
            and (r1_parts[-2] == "fq" or r1_parts[-2] == "fastq")
        )
        else None
    )
    if prefix_pos == None:
        prefix_pos = -1 if ((r1_parts[-1] == "fq" or r1_parts[-1] == "fastq")) else None
    return prefix_pos


def run_trim(read_list, output_dir, job_data, tool_params):
    rcount = 0
    trimmed_reads = []
    threads = get_param("trim_galore", "-p", tool_params)
    trim_base_cmd = ["trim_galore", "-j", str(threads) if threads else "1"]
    for r in read_list:
        tr = copy.deepcopy(r)
        tr["fastqc"] = []
        cur_check = []
        rename_files = {}
        rcount += 1
        # trim_galore --gzip --paired -o ../fonsynbiothr/trimmed_reads/ ../fonsynbiothr/fastq_files/926M_RNA_S8_L001_R1_001.fastq ../fonsynbiothr/fastq_files/926M_RNA_S8_L001_R2_001.fastq
        trim_cmd = trim_base_cmd + ["--gzip", "-o", output_dir]
        if "read2" in r:
            trim_cmd += ["--paired", r["read1"], r["read2"]]
            pre_pos = find_prefix(r["read1"])
            old_name = os.path.join(
                output_dir,
                ".".join(os.path.basename(r["read1"]).split(".")[0:pre_pos])
                + "_val_1.fq.gz",
            )
            cur_check.append(old_name)
            new_name = os.path.join(
                output_dir,
                ".".join(os.path.basename(r["read1"]).split(".")[0:pre_pos])
                + "_ptrim.fq.gz",
            )
            rename_files[old_name] = new_name
            tr["read1"] = new_name
            pre_pos = find_prefix(r["read2"])
            old_name = os.path.join(
                output_dir,
                ".".join(os.path.basename(r["read2"]).split(".")[0:pre_pos])
                + "_val_2.fq.gz",
            )
            cur_check.append(old_name)
            new_name = os.path.join(
                output_dir,
                ".".join(os.path.basename(r["read2"]).split(".")[0:pre_pos])
                + "_ptrim.fq.gz",
            )
            rename_files[old_name] = new_name
            tr["read2"] = new_name
        else:
            trim_cmd += [r["read1"]]
            pre_pos = find_prefix(r["read1"])
            old_name = os.path.join(
                output_dir,
                ".".join(os.path.basename(r["read1"]).split(".")[0:pre_pos])
                + "_trimmed.fq.gz",
            )
            cur_check.append(old_name)
            new_name = os.path.join(
                output_dir,
                ".".join(os.path.basename(r["read1"]).split(".")[0:pre_pos])
                + "_strim.fq.gz",
            )
            rename_files[old_name] = new_name
            tr["read1"] = new_name
        print((" ".join(trim_cmd)))
        subprocess.run(trim_cmd, check=True)
        check_passed = True
        for c in cur_check:
            check_passed = check_passed and os.path.exists(c)
            if c in rename_files:
                os.rename(c, rename_files[c])
        if check_passed:
            trimmed_reads.append(tr)
        else:
            sys.stderr.write("Trimming reads failed at " + " ".join(cur_check))
            sys.exit()
    return trimmed_reads

def get_minimap_preset(platform=None):
    """
    Try to determine the preset string, based on the
    jobspec reads platform, if present.
    default: []
    """
    lookup_dict = {
        "pacbio":            "map-pb",    # PacBio CLR genomic reads
        "nanopore":          "map-ont",   # Oxford Nanopore genomic reads
        "pacbio_hifi":       "map-hifi",  # PacBio HiFi/CCS genomic reads (v2.19 or later)
        "packbio_hifi_2_18": "asm20",     # PacBio HiFi/CCS genomic reads (v2.18 or earlier)
        "illumina":          "sr",        # short genomic (possibly paired-end) reads
    }
    if platform in lookup_dict:
        return ["-x", lookup_dict[platform]]
    else:
        return []

def run_alignment(genome_list, read_list, parameters, output_dir, job_data):
    # modifies condition_dict sub replicates to include 'bowtie' dict recording output files

    # Get a minimap2 command ready for potential duty.
    mm2_cmd = ["minimap2", "-a"] # . . . '-a' is to output in SAM format
    mm2_params = parameters.get("minimap2", {})
    # thread count (tool default is 3)
    if "-t" in mm2_params:
        mm2_cmd += ["-t", str(mm2_params["-t"])]

    for genome in genome_list:
        genome_link = genome["genome_link"]
        final_cleanup = []
        hisat2_used = False
        if "hisat_index" in genome and genome["hisat_index"]:
            archive = tarfile.open(genome["hisat_index"])
            archive.extractall(path=output_dir)
            indices = []
            for ht2_filename in archive.getnames():
                home_name = os.path.join(output_dir, os.path.basename(ht2_filename))
                if not os.path.exists(home_name):
                    shutil.move(os.path.join(output_dir, ht2_filename), home_name)
                indices.append(home_name)
            final_cleanup += indices
            index_prefix = os.path.join(
                output_dir, re.sub(r"\.?[0-9]*\.ht2$", "", os.path.basename(indices[0]))
            )
            archive.close()
            # index_prefix = os.path.join(
            #     output_dir,
            #     os.path.basename(genome["hisat_index"]).replace(".ht2.tar", ""),
            # )  # somewhat fragile convention. tar prefix is underlying index prefix
            cmd = ["hisat2", "--dta-cufflinks", "-x", index_prefix]
            hisat2_params = parameters.get("hisat2", {})
            # thread count
            if "-p" in hisat2_params:
                cmd += ["-p", str(hisat2_params["-p"])]
            hisat2_used = True
        else:
            subprocess.run(["bowtie2-build", genome_link, genome_link], check=True)
            bt2_threads = parameters.get("bowtie2", {}).get("-p", BT2_THREADS)
            cmd = ["bowtie2", "-p", str(bt2_threads), "-x", genome_link]

        print(f"{cmd=}")

        target_dir = genome["output"]
        for r in read_list:
            rcount = 0
            cur_cleanup = []
            rcount += 1
            samstat_cmd = ["samstat"]
            cur_cmd = list(cmd)  # hisat2 or bowtie2
            read2 = False
            # generate output file name, based on input read file(s) name(s).
            out_name = "_".join([Path(i[1]).stem.replace(" ", "") for i in r.items() if "read" in i[0]])
            sam_path = Path(target_dir) / f"{out_name}.sam"
            unmapped_fq_gz_path = Path(target_dir) / f"{out_name}.unmapped.fq.gz"
            if "read2" in r:  # illumina, using hisat2 or bowtie2
                cur_cmd += [
                    "-1", link_space(r["read1"]), "-2", link_space(r["read2"]),
                    "-S", str(sam_path),
                    "--un-conc-gz", str(unmapped_fq_gz_path),
                ]
                read2 = True
            else:
                if hisat2_used:
                    cur_cmd += [
                        "-U", link_space(r["read1"]),
                        "-S", str(sam_path),
                        "--un-gz", str(unmapped_fq_gz_path),
                    ]
                else:  # minimap2
                    cur_cmd = mm2_cmd
                    cur_cmd += [
                        *get_minimap_preset(r.get("platform", None)),
                        "-o", str(sam_path),
                        genome_link,
                        link_space(r["read1"]),
                        ]
            print(f"{cur_cmd=}")
            cur_cleanup.append(str(sam_path))
            bam_file_all =        str(sam_path.with_suffix(".all.bam"))
            bam_file_aligned =    str(sam_path.with_suffix(".aligned.bam"))
            bam_file_sort_name =  str(sam_path.with_suffix(".aligned.name.bam"))
            fastq_file_aligned =  str(sam_path.with_suffix(".aligned.fq"))
            fastq_file_aligned2 = str(sam_path.with_suffix(".aligned.2.fq"))
            samstat_cmd.append(bam_file_all)
            # make bam file information available for subsequent steps
            r[genome["genome"]] = {}
            r[genome["genome"]]["bam"] = bam_file_aligned
            if os.path.exists(bam_file_aligned):
                sys.stderr.write(
                    bam_file_aligned + " alignments file already exists. skipping\n"
                )
            else:
                subprocess.run(cur_cmd, check=True)  # call mapping tool
                view_threads = parameters.get("samtools_view", {}).get("-p", SAM_THREADS)
                with open(bam_file_all, "w") as outstream:
                    samtools_view = subprocess.Popen(
                        ["samtools", "view", "-@", str(view_threads), "-Su", str(sam_path)],
                        stdout=subprocess.PIPE,
                    )
                    samtools_sort = subprocess.Popen(
                        ["samtools", "sort", "-o", "-", "-"],
                        stdin=samtools_view.stdout,
                        stdout=outstream,
                    )
                    samtools_sort.communicate()

                # minimap2 can't write the unmapped reads like bowtie2 and hisat2
                if not (read2 or hisat2_used): # minimap2
                    with unmapped_fq_gz_path.open('w') as un_fq_gz_hdl:
                        subprocess.run(
                            [
                                "samtools", "fastq",
                                "--threads", str(view_threads),
                                "-n", "-f", "4",
                                bam_file_all
                            ],
                            check=True,
                            stdout=un_fq_gz_hdl,
                        )

                with open(bam_file_aligned, "w") as outstream:
                    subprocess.run(
                        [
                            "samtools",
                            "view",
                            "-@",
                            str(view_threads),
                            "-b",
                            "-F",
                            "4",
                            bam_file_all,
                        ],
                        check=True,
                        stdout=outstream,
                    )
                sam_sort = parameters.get("samtools_sort", {}).get("-p", "1")
                if not sam_sort:
                    sam_sort = SAM_THREADS
                sort_name_cmd = [
                    "samtools",
                    "sort",
                    "-n",
                    "-@",
                    str(sam_sort),
                    bam_file_aligned,
                ]
                bam2fq_cmd = [
                    "bedtools",
                    "bamtofastq",
                    "-i",
                    bam_file_sort_name,
                    "-fq",
                    fastq_file_aligned,
                ]
                bam2fqgz_cmd = ["gzip", fastq_file_aligned]
                if read2:  # paired end
                    bam2fq_cmd += ["-fq2", fastq_file_aligned2]
                    bam2fqgz_cmd += [fastq_file_aligned2]
                print(" ".join(sort_name_cmd))
                # Need to sort by name to convert to fastq: samtools sort -n myBamFile.bam myBamFile.sortedByName
                subprocess.run(
                    sort_name_cmd, check=True, stdout=open(bam_file_sort_name, "w")
                )
                cur_cleanup.append(bam_file_sort_name)
                print((" ".join(bam2fq_cmd)))
                with open(os.path.join(target_dir, "bedtools.log.txt"), "a") as fd:
                    print("Redirecting bedtools stderr to a log.", file=sys.stderr)
                    print(bam_file_sort_name, file=fd)
                    fd.flush()
                    subprocess.run(bam2fq_cmd, check=True, stderr=fd)
                subprocess.run(bam2fqgz_cmd, check=True)
                samtools_index_threads = parameters.get("samtools_index", {}).get(
                    "-p", "1"
                )
                if not samtools_index_threads:
                    samtools_index_threads = SAM_THREADS
                subprocess.run(
                    [
                        "samtools",
                        "index",
                        "-@",
                        str(samtools_index_threads),
                        bam_file_aligned,
                    ],
                    check=True,
                )
                print((" ".join(samstat_cmd)))
                try:
                    subprocess.run(samstat_cmd, check=True)
                except subprocess.CalledProcessError as cpe:
                    with open(bam_file_all + ".samstat.html", "w") as samstat_html:
                        samstat_html.write("<!DOCTYPE html>\n<html>\n<body>\n")
                        samstat_html.write(
                            "<h2>Your samstat job had a segfault error. The results could not be displayed.</h2>\n"
                        )
                        samstat_html.write(
                            "<p><b>return code:</b> {}</p>\n<p><b>command:</b> {}</p>\n<p><b>output:</b> {}</p>\n<p><b>stdout:</b> {}</p>\n<p><b>stderr:</b> {}</p>\n".format(
                                cpe.returncode,
                                cpe.cmd,
                                cpe.output,
                                cpe.stdout,
                                cpe.stderr,
                            ),
                        )
                        samstat_html.write("</html>\n</body>\n")
                cur_cleanup.append(bam_file_all)
            for garbage in cur_cleanup:
                if Path(garbage).exists(): os.remove(garbage)
        files = [
            f
            for f in os.listdir(target_dir)
            if os.path.isfile(os.path.join(target_dir, f))
        ]
        for f in files:
            if f.endswith(".fq.1.gz"):
                os.rename(
                    os.path.join(target_dir, f),
                    os.path.join(target_dir, f[:-8] + ".1.fq.gz"),
                )
            elif f.endswith(".fq.2.gz"):
                os.rename(
                    os.path.join(target_dir, f),
                    os.path.join(target_dir, f[:-8] + ".2.fq.gz"),
                )
            elif f.endswith(".fq.1"):
                os.rename(
                    os.path.join(target_dir, f),
                    os.path.join(target_dir, f[:-5] + ".1.fq"),
                )
            elif f.endswith(".fq.2"):
                os.rename(
                    os.path.join(target_dir, f),
                    os.path.join(target_dir, f[:-5] + ".2.fq"),
                )
        for garbage in final_cleanup:
            os.remove(garbage)


def unzip(path, value):
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as f_in:
            with open(path[0 : len(path) - 3], "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        # subprocess.run(["gunzip", path])
        value.value = 3


def paired_filter(read_list, parameters, output_dir, job_data):
    print("Running paired filter.")
    for r in read_list:
        if "read2" in r:
            value1 = multiprocessing.Value("i", 0, lock=False)
            value2 = multiprocessing.Value("i", 0, lock=False)
            p1 = Process(target=unzip, args=(r["read1"], value1))
            p2 = Process(target=unzip, args=(r["read2"], value2))
            p1.start()
            p2.start()
            p1.join()
            p2.join()
            r["read1"] = r["read1"][0 : len(r["read1"]) - value1.value]
            r["read2"] = r["read2"][0 : len(r["read2"]) - value2.value]
            pair_cmd = ["fastq_pair", r["read1"], r["read2"]]
            subprocess.run(pair_cmd)
            r["read1"] += ".paired.fq"
            r["read2"] += ".paired.fq"
    return read_list

def run_hostile(read_list, output_dir, job_data, tool_params):
    """
    Run Hostile to remove host sequences from short and long read (meta)genomes.
    https://github.com/bede/hostile
    """
    output_dir = Path(output_dir)
    print(f"{output_dir=}", file=sys.stderr)

    # Set the cache dir for the index files
    # See https://github.com/bede/hostile/issues/32
    cache_dir = output_dir / "hostile"
    cache_dir.mkdir(exist_ok=True, parents=True)
    os.environ["HOSTILE_CACHE_DIR"] = str(cache_dir)
    print(f"{os.environ['HOSTILE_CACHE_DIR']=}", file=sys.stderr)

    log_path = output_dir / "hostile_report.txt"
    print(f"{log_path=}", file=sys.stderr)

    base_cmd = ["hostile", "clean", "--force",
           "--out-dir", str(output_dir), "--fastq1"]

    for r in read_list:
        print(f"{r=}", file=sys.stderr)
        if "read2" in r:
            cmd = base_cmd + [r["read1"], "--fastq2", r["read2"]]
        else:
            cmd = base_cmd + [r["read1"],]

        with log_path.open('w') as log_hdl:
            subprocess.run(cmd, check=True, stdout=log_hdl)

def get_genome(parameters, host_manifest={}):
    target_file = os.path.join(parameters["output_path"], parameters["gid"] + ".fna")
    genome = {"genome_link": target_file}
    print("GID: {}".format(parameters["gid"]), file=sys.stderr)
    if not os.path.exists(target_file):
        taxid = str(parameters["gid"]).split(".")[0]
        if taxid in host_manifest:
            # Only the index is needed.
            # genome_url = host_manifest[taxid]["patric_ftp"] + "_genomic.fna"
            # print(genome_url)
            # with closing(request.urlopen(genome_url)) as r:
            #     with open(target_file, "wb") as f:
            #         shutil.copyfileobj(r, f)
            index_url = host_manifest[taxid]["patric_ftp"] + "_genomic.ht2.tar"
            index_file = os.path.join(
                parameters["output_path"], parameters["gid"] + ".ht2.tar"
            )
            print(index_url)
            with closing(request.urlopen(index_url)) as r:
                with open(index_file, "wb") as f:
                    shutil.copyfileobj(r, f)
            genome["hisat_index"] = index_file
        else:
            # parameters["data_api"], parameters["gid"])
            genome_url = (
                "{data_api}/genome_sequence/?eq(genome_id,{gid})&limit(25000)".format(
                    **parameters
                )
            )
            headers = {"accept": "application/sralign+dna+fasta"}
            req = requests.Request("GET", genome_url, headers=headers)
            print(genome_url)
            # print "switch THE HEADER BACK!"
            # headers = {'Content-Type': 'application/x-www-form-urlencoded; charset=utf-8'}
            authenticateByEnv(req)
            prepared = req.prepare()
            # pretty_print_POST(prepared)
            s = requests.Session()
            response = s.send(prepared)
            handle = open(target_file, "wb")
            if not response.ok:
                sys.stderr.write("API not responding. Please try again later.\n")
                sys.exit(2)
            else:
                for block in response.iter_content(1024):
                    handle.write(block)
    sys.stdout.flush()
    return genome

sra_to_local_pltfrm = {
    "ILLUMINA": "illumina",
    "OXFORD_NANOPORE": "nanopore",
    "PACBIO_SMRT": "pacbio",
    "ION_TORRENT": None,
    "SOLID": None,
    "LS454": None,
    "CAPILLARY": None,
    }

def setup(job_data, output_dir, tool_params):
    genome_ids = []
    ref_id = job_data.get("reference_genome_id", None)
    if ref_id != None:
        genome_ids.append(ref_id)
    if genome_ids:
        manifest = getHostManifest()
    genome_list = []
    for gid in genome_ids:
        genome = {}
        # cheat. this will need expansion if you want to support multiple genomes
        job_data["gid"] = gid
        genome.update(get_genome(job_data, manifest))
        genome["gid"] = gid
        genome["genome"] = gid
        genome["output"] = output_dir
        genome_list.append(genome)
    read_list = []
    rcount = 0
    for r in (
        job_data.get("paired_end_libs", [])
        + job_data.get("single_end_libs", [])
        + job_data.get("srr_libs", [])
    ):
        if "read" in r:
            r["read1"] = r.pop("read")
        read_list.append(r)
        r["fastqc"] = []
        target_dir = output_dir
        rcount += 1
        if "srr_accession" in r:
            srr_id = r["srr_accession"]
            meta_file = os.path.join(target_dir, srr_id + "_meta.txt")
            sra_cmd = [
                "p3-sra",
                "--out",
                target_dir,
                "--metadata-file",
                meta_file,
                "--id",
                srr_id,
            ]
            print((" ".join(sra_cmd)))
            subprocess.run(sra_cmd, check=True)
            with open(meta_file) as f:
                job_meta = json.load(f)
                if "platform_name" in job_meta[0]:
                    r.setdefault("platform", sra_to_local_pltfrm[job_meta[0]["platform_name"]])
                files = job_meta[0].get("files", [])
                if len(files) > 0:
                    for _, f in enumerate(files):
                        if f.endswith("_2.fastq"):
                            r["read2"] = os.path.join(target_dir, f)
                        elif f.endswith("_1.fastq"):
                            r["read1"] = os.path.join(target_dir, f)
                        elif f.endswith(".fastq"):
                            r["read1"] = os.path.join(target_dir, f)
                        elif f.endswith("fastqc.html"):
                            r["fastqc"].append(os.path.join(target_dir, f))
    recipe = job_data.get("recipe", [])
    new_read_list = []
    for r in read_list:
        if "read1" in r and "read2" not in r and os.path.exists(r["read1"]):
            r["read1"] = moveRead(r["read1"])
            new_read_list.append(r)
        if (
            "read1" in r
            and "read2" in r
            and os.path.exists(r["read1"])
            and os.path.exists(r["read2"])
        ):
            r["read1"] = moveRead(r["read1"])
            r["read2"] = moveRead(r["read2"])
            new_read_list.append(r)
    return genome_list, new_read_list, recipe


def gzipMove(source, dest):
    with open(source, "rb") as f_in:
        with gzip.open(dest, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def moveRead(filepath):
    directory, base = os.path.split(filepath)
    # Don't try to escape to prevent escape shenanigans.
    new_base = re.sub(r"[^a-zA-Z0-9_\-\.]", "_", base)
    # new_base = re.sub(r"[\`~\"'!@#$%^&*(){}\[\]|<>]", "_", base)
    # new_base = base
    # for spec_char, repl in SPECIAL_CHARS:
    #     new_base = new_base.replace(spec_char, repl)
    if new_base != base:
        newpath = os.path.join(directory, new_base)
        os.rename(filepath, newpath)
        return newpath
    else:
        return filepath

def run_fq_util(job_data, output_dir, tool_params={}):
    # arguments:
    # list of genomes [{"genome":somefile,"annotation":somefile}]
    # dictionary of library dictionaries structured as {libraryname:{library:libraryname, replicates:[{read1:read1file, read2:read2file}]}}
    # parametrs_file is a json keyed parameters list.
    # Example tool_params: '{"fastqc":{"-p":"2"},"trim_galore":{"-p":"2"},"bowtie2":{"-p":"2"},"hisat2":{"-p":"2"},"samtools_view":{"-p":"2"},"samtools_index":{"-p":"2"}}'
    output_dir = os.path.abspath(output_dir)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    genome_list, read_list, recipe = setup(job_data, output_dir, tool_params)
    # print("genome_list: {}\nread_list: {}\nrecipe: {}".format(genome_list, read_list, recipe), file=sys.stdout)
    # sys.stdout.flush()
    for step in recipe:
        step = step.upper()
        if step == "TRIM":
            trimmed_reads = run_trim(read_list, output_dir, job_data, tool_params)
            read_list = trimmed_reads
        elif step == "PAIRED_FILTER":
            read_list = paired_filter(read_list, tool_params, output_dir, job_data)
        elif step == "FASTQC":
            run_fastqc(read_list, output_dir, job_data, tool_params)
        elif step == "ALIGN":
            run_alignment(genome_list, read_list, tool_params, output_dir, job_data)
        elif step == "SCRUB_HUMAN":
            run_hostile(read_list, output_dir, job_data, tool_params)
        else:
            print("Skipping step. Not found: {}".format(step), file=sys.stderr)
    if len(recipe) == 1 and recipe[0].upper() == "PAIRED_FILTER":
        for r in read_list:
            if "read2" in r:
                dest1 = os.path.join(output_dir, os.path.basename(r["read1"]) + ".gz")
                dest2 = os.path.join(output_dir, os.path.basename(r["read2"]) + ".gz")
                p1 = Process(target=gzipMove, args=(r["read1"], dest1))
                p2 = Process(target=gzipMove, args=(r["read2"], dest2))
                p1.start()
                p2.start()
                p1.join()
                p2.join()
