"""
python module with various methods for bacterial annotation and comparative analysis
"""

from __future__ import print_function
import sys,os,subprocess,glob,shutil
from Bio import Entrez
Entrez.email = 'A.N.Other@example.com'
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Phylo, AlignIO
import matplotlib as mpl
import pylab as plt
import seaborn as sns
from matplotlib.colors import ListedColormap, LogNorm
import numpy as np
import pandas as pd

def fastq_to_dataframe(f, size=None):
    """Convert fastq to dataframe.
        size: limit to the first reads of total size
        Returns: dataframe with reads
    """

    import HTSeq
    ext = os.path.splitext(f)[1]
    if ext=='.fastq' or ext=='.gz':
        ffile = HTSeq.FastqReader(f, "solexa")
    elif ext == '.fa':
        ffile = HTSeq.FastaReader(f)
    else:
        return
    if size != None:
        sequences = [(s.name, s.seq, s.descr) for s in islice(fastfile, i, i+size)]
    else:
        sequences = [(s.name,s.seq) for s in ffile]
    df = pd.DataFrame(sequences,columns=['id','seq'])
    return df

def dataframe_to_fasta(df, seqkey='translation', idkey='locus_tag',
                     descrkey='description',
                     outfile='out.faa'):
    """Genbank features to fasta file"""

    seqs=[]
    for i,row in df.iterrows():
        if descrkey in df.columns:
            d=row[descrkey]
        else:
            d=''
        rec = SeqRecord(Seq(row[seqkey]),id=row[idkey],
                            description=d)
        seqs.append(rec)
    SeqIO.write(seqs, outfile, "fasta")
    return outfile

def fasta_to_dataframe(infile, header_sep=None, key='name', seqkey='sequence'):
    """Get fasta proteins into dataframe"""

    recs = SeqIO.parse(infile,'fasta')
    keys = [key,seqkey,'description']
    data = [(r.name,str(r.seq),str(r.description)) for r in recs]
    df = pd.DataFrame(data,columns=(keys))
    df['type'] = 'CDS'
    #fix bad names
    if header_sep not in ['',None]:
        df[key] = df[key].apply(lambda x: x.split(header_sep)[0],1)
    df[key] = df[key].str.replace('|','_')
    return df

def get_blast_results(filename):
    """
    Get blast results into dataframe. Assumes column names from local_blast method.
    Returns:
        dataframe
    """

    cols = ['qseqid','sseqid','qseq','sseq','pident','qcovs','length','mismatch','gapopen',
            'qstart','qend','sstart','send','evalue','bitscore','stitle']
    res = pd.read_csv(filename, names=cols, sep='\t')
    #res = res[res['pident']>=ident]
    return res

def local_blast(database, query, output=None, maxseqs=50, evalue=0.001,
                    compress=False, cmd='blastn', cpus=4, show_cmd=False, **kwargs):
    """Blast a local database.
    Args:
        database: local blast db name
        query: sequences to query, list of strings or Bio.SeqRecords
    Returns:
        pandas dataframe with top blast results
    """

    if output == None:
        output = os.path.splitext(query)[0]+'_blast.txt'
    from Bio.Blast.Applications import NcbiblastxCommandline
    outfmt = '"6 qseqid sseqid qseq sseq pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore stitle"'
    cline = NcbiblastxCommandline(query=query, cmd=cmd, db=database,
                                 max_target_seqs=maxseqs,
                                 outfmt=outfmt, out=output,
                                 evalue=evalue, num_threads=cpus, **kwargs)
    if show_cmd == True:
        print (cline)
    stdout, stderr = cline()
    return

def blast_sequences(database, seqs, labels=None, **kwargs):
    """
    Blast a set of sequences to a local or remote blast database
    Args:
        database: local or remote blast db name
                  'nr', 'refseq_protein', 'pdb', 'swissprot' are valide remote dbs
        seqs: sequences to query, list of strings or Bio.SeqRecords
        labels: list of id names for sequences, optional but recommended
    Returns:
        pandas dataframe with top blast results
    """

    remotedbs = ['nr','refseq_protein','pdb','swissprot']
    res = []
    if not type(seqs) is list:
        seqs = [seqs]
    if labels is None:
        labels = seqs
    recs=[]
    #print (labels)
    for seq, name in zip(seqs,labels):
        if type(seq) is not SeqRecord:
            rec = SeqRecord(Seq(seq),id=name)
        else:
            rec = seq
            name = seq.id
        recs.append(rec)
    SeqIO.write(recs, 'tempseq.fa', "fasta")
    if database in remotedbs:
        remote_blast(database, 'tempseq.fa', **kwargs)
    else:
        local_blast(database, 'tempseq.fa', **kwargs)
    df = get_blast_results(filename='tempseq_blast.txt')
    return df

def clustal_alignment(filename=None, seqs=None, command="clustalw"):
    """Align 2 sequences with clustal"""

    if filename == None:
        filename = 'temp.faa'
        SeqIO.write(seqs, filename, "fasta")
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import ClustalwCommandline
    cline = ClustalwCommandline(command, infile=filename)
    stdout, stderr = cline()
    align = AlignIO.read(name+'.aln', 'clustal')
    return align

def muscle_alignment(filename=None, seqs=None):
    """Align 2 sequences with muscle"""

    if filename == None:
        filename = 'temp.faa'
        SeqIO.write(seqs, filename, "fasta")
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input=filename, out=name+'.txt')
    stdout, stderr = cline()
    align = AlignIO.read(name+'.txt', 'fasta')
    return align

def align_nucmer(file1, file2):
    cmd='nucmer --maxgap=500 --mincluster=100 --coords -p nucmer %s %s' %(file1, file2)
    print (cmd)
    subprocess.check_output(cmd,shell=True)
    df = read_nucmer_coords('nucmer.coords')
    return df

def read_nucmer_coords(cfile):
    cols=['S1','E1','S2','E2','LEN 1','LEN 2','IDENT','TAG1','TAG2']
    a=pd.read_csv(cfile,sep='[\s|]+',skiprows=5,names=cols,engine='python')
    a = a.sort_values(by='TAG2',ascending=False)
    return a

def align_reads(file1, file2, idx, out):
    """align reads to ref"""

    cmd = 'bwa mem -M -t 8 %s %s %s | samtools view -bS - > %s' %(idx,files[0],files[1],out)
    if not os.path.exists(out):
        print (cmd )
        subprocess.check_output(cmd, shell=True)
    return

def align_info(bamfile):
    cmd = 'samtools flagstat %s' %bamfile
    temp=subprocess.check_output(cmd, shell=True)
    print (temp)
    return

def variants_call(name, ref, out):
    bamfile = '%s/%s.bam' %(out,name)
    cmd = 'samtools sort {b} > {b}.sorted && samtools index {b}.sorted'.format(b=bamfile)
    print (cmd)
    #subprocess.check_output(cmd, shell=True)
    cmd = 'samtools mpileup -uf genomes/{r}.fa {b}.sorted | bcftools call -mv \
    > {o}/{n}.vcf'.format(b=bamfile,n=name,r=ref,o=out)
    print (cmd)
    #subprocess.check_output(cmd, shell=True)
    cmd = 'bedtools intersect -a {gff} -b {o}/{n}.vcf -wa -u > {o}/{n}_variants.bed'.format(n=name,r=ref,gff=gff,o=out)
    print (cmd)

def search_genbank(term='', filt=None):
    request = Entrez.esearch(db="nuccore", term=term, field="title", FILT=filt, rettype='xml')
    result = Entrez.read(request)
    idlist = result['IdList']
    return idlist

#phylo trees

def get_tree(aln, kind='nj'):
    from Bio.Phylo.TreeConstruction import DistanceCalculator,DistanceTreeConstructor
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    return dm, tree

def clear_clades(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if 'Inner' in clade.name :
            clade.name = ''

        names[clade.name] = clade
    return names

def draw_tree(tree, root=None, labels=None, clear=True, title=''):
    f,ax=plt.subplots(figsize=(10,8))
    if clear == True:
        try:
            clear_clades(tree)
        except:
            pass
    if root != None:
        tree.root_with_outgroup(root)
    if labels != None:
        for clade in tree.get_terminals():
            key = clade.name
            if key in labels:
                clade.name = '%s; %s' %(key,labels[key])
                #clade.name = labels[key]

    Phylo.draw(tree,axes=ax,axis=('off',),branch_labels=None,show_confidence=False)
    #Phylo.draw_graphviz(tree,axes=ax,node_size=0)
    ax.set_title(title,fontsize=16)
    return f, tree

def show_alignment(aln, diff=False, offset=0):
    """
    Show a sequence alignment
        Args:
            aln: alignment
            diff: whether to show differences
    """

    ref = aln[0]
    l = len(aln[0])
    n=60
    chunks = [(i,i+n) for i in range(0, l, n)]
    for c in chunks:
        start,end = c
        lbls = np.arange(start,end,10)-offset
        print (('%-21s' %'name'),''.join([('%-10s' %i) for i in lbls]))
        print (('%21s' %ref.id[:20]), ref.seq[start:end])

        if diff == True:
            for a in aln[1:]:
                diff=''
                for i,j in zip(ref,a):
                    if i != j:
                        diff+=j
                    else:
                        diff+='-'
                name = a.id[:20]
                print (('%21s' %name), diff[start:end])
        else:
            for a in aln[1:]:
                name = a.id[:20]
                print (('%21s' %name), a.seq[start:end])
    return

def get_assembly_summary(id):

    from Bio import Entrez
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record

def get_assemblies(term=None, ids=None, download=True, path='assemblies'):
    """Download genbank assemblies for a given search term.
    Args:
        ids: entrez ids for assemblies
        term: usually organism name"""

    from Bio import Entrez
    Entrez.email = "A.N.Other@example.com"
    if ids == None:
        handle = Entrez.esearch(db="assembly", term=term, retmax='200')
        record = Entrez.read(handle)
        ids = record['IdList']
        print (f'found {len(ids)} ids')
    links = []   
    for id in ids:
        #get summary
        summary = get_assembly_summary(id)
        #get ftp link
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            continue
        label = os.path.basename(url)
        link = os.path.join(url,label+'_genomic.fna.gz')
        print (link)
        links.append(link)
        if download == True:
            #download link
            urllib.request.urlretrieve(link, f'{label}.fna.gz')
    return links

def get_url_from_path(url):
    """get full path for genomic fasta from url"""
    
    label = os.path.basename(url)
    link = os.path.join(url,label+'_genomic.fna.gz')
    return link

 