import re

class Style_default:
    ID = 'ID'
    ATTRIBUTES = ['ID', 'Parent']
    ATTRIBUTE_DELIMITERS = ','

    @staticmethod
    def is_gene(feature):
        return feature.feature == 'gene'

class Style_phytozome(Style_default):
    ID = 'Name'

class Style_Refseq(Style_default):
    ATTRIBUTES = ['ID', 'Parent', 'gene_biotype']
    ATTRIBUTE_DELIMITERS = r',|\|'

    @staticmethod
    def is_gene(feature):
        return feature.feature == 'gene' and 'gene_biotype' in feature.attrs and feature.attrs['gene_biotype'] == ['protein_coding'] 


class Feature:
    def __init__(self, raw, style, idGen=None, saveRAM=True):
        self.raw = raw
        self.style = style
        # seqid - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seq ID must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
        # source - name of the program that generated this feature, or the data source (database or project name)
        # feature - type of feature. Must be a term or accession from the SOFA sequence ontology
        # start - Start position of the feature, with sequence numbering starting at 1.
        # end - End position of the feature, with sequence numbering starting at 1.
        # score - A floating point value.
        # strand - defined as + (forward) or - (reverse).
        # phase - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
        # attributes - A semicolon-separated list of tag-value pairs, providing additional information about each feature. Some of these tags are predefined, e.g. ID, Name, Alias, Parent - see the GFF documentation for more details.
        self.seqid, self.source, self.feature, self.start, self.end, self.score, self.strand, self.phase, self.attributes = tuple(raw.split('\t'))
        self.start = int(self.start)
        self.end = int(self.end)
        self.size = self.end - self.start + 1
        self.attrs = {
            a[0]: re.split(self.style.ATTRIBUTE_DELIMITERS, a[1]) 
                for a in [attribute.split('=') 
                    for attribute in self.attributes.split(';') if '=' in attribute]
                if a[0] in style.ATTRIBUTES
            }
        self.id = self.attrs[style.ID][0] if style.ID in self.attrs else idGen()
        self.name = self.attrs['Name'][0] if 'Name' in self.attrs else None
        self.parents = self.attrs['Parent'] if 'Parent' in self.attrs else []
        self.parsed_parents = {}
        if saveRAM:
            del self.attributes
    
    def toStore(self):
        return '\t'.join(
            [self.seqid, self.source, self.feature, str(self.start), str(self.end), self.score, self.strand, self.phase, 
            'ID={};{}'.format(self.id, ('Parent=' + ','.join(self.parents) + ';') if len(self.parents) > 0 else '')]
        )
    
    def parsed_attributes(self):
        if not self.attributes:
            raise Exception('Has not attributes, loaded with save RAM?')
        return {
            a[0]: re.split(self.style.ATTRIBUTE_DELIMITERS, a[1]) 
                for a in [attribute.split('=') 
                    for attribute in self.attributes.split(';') if '=' in attribute]
            }
    
    def is_gene(self):
        return self.style.is_gene(self)
    
    def __str__(self):
        return "[{} {}] {} {}({}:{}){}".format(
            self.feature, self.id, self.seqid, self.strand, self.start, self.end, ' <= ' + ', '.join(self.parents) if self.parents else ''
            )
    
    def __repr__(self):
        return "{}: {} ({})\n{}: {}{}:{}\n{}".format(
            self.feature, self.id, ' <= ' + ', '.join(self.parents) if self.parents else '', 
            self.seqid, self.strand, self.start, self.end, '; '.join(['{}: {}'.format(k, ', '.join(v)) for k, v in self.attrs.items() if not k in [self.style.ID, 'Parent']]) 
        )

class Gene(Feature):
    def __init__(self, raw, style, idGen):
        Feature.__init__(self, raw, style, idGen=idGen)
        self.mrnas = []
    
    def toStore(self):
        return Feature.toStore(self) + (('\n' + '\n'.join([mrna.toStore() for mrna in self.mrnas])) if len(self.mrnas) > 0 else '')

class Mrna(Feature):
    def __init__(self, raw, gene, style, idGen):
        Feature.__init__(self, raw, style, idGen=idGen)
        self.gene = gene
        gene.mrnas.append(self)
        self.exons = []
        self.introns = []
        self.five_prime_UTRs = []
        self.three_prime_UTRs = []
        self.cds = []

    def generateIntrons(self):
        if len(self.exons) > 1 and len(self.introns) < 1:
            exons = list(sorted(self.exons, key=lambda e: e.start))
            for e1, e2 in [(exons[i], exons[i+1]) for i in range(len(exons)-1)]:
                self.introns.append(Feature(
                    raw='\t'.join([e1.seqid, 'gen', 'intron',str( e1.end + 1), str(e2.start - 1), '.', e1.strand, '.', 
                    'ID={};Parent={}'.format(e1.id +'-'+e2.id, self.id)]), 
                    style=e1.style
                    ))
        return self.introns

    def toStore(self):
        features = self.exons + self.introns + self.five_prime_UTRs + self.three_prime_UTRs + self.cds
        return Feature.toStore(self) + (('\n' + '\n'.join([
            feature.toStore() for feature in features
            ])) if len(features) > 0 else '')


class GFF:

    def debug_info(self, txt):
        if 'I' in self.debug:
            print('[INFO] ' + txt)

    def debug_warning(self, txt):
        if 'W' in self.debug:
            print('[WARN] ' + txt)

    def debug_error(self, txt):
        if 'E' in self.debug:
            print('[ERRO] ' + txt)

    def idGenerator(self):
        self.gen += 1
        return 'ID_gen_%d' % self.gen

    def __init__(self, file, target_genes=None, target_seqs=None, style=Style_default, debug=['E'], saveRAM=True):
        self.file_name = file
        self.debug = debug
        self.style = style
        self.gen = 0

        if not type(self.style) in [type(s) for s in [Style_default, Style_Refseq]]:
            self.debug_error('Style must be provided.')
            raise Exception('Style must be one of Style_default or Style_Refseq')
        
        self.target_genes = [target_genes] if type(target_genes) == str else target_genes
        self.target_seqs = [target_seqs] if type(target_seqs) == str else target_seqs

        if not self.target_genes is None:
            self.debug_info('Filtering for genes: ' + ', '.join(self.target_genes))
        
        if not self.target_seqs is None:
            self.debug_info('Filtering for seqs: ' + ', '.join(self.target_seqs))

        if type(self.target_seqs) == type(self.target_genes) == list and len(self.target_seqs) > 0 and len(self.target_genes) > 0:
            self.debug_warning('Ambigous filter definition privilegies inclusion.')

        self.debug_info('Opening file {}'.format(self.file_name))

        try:
            file = open(self.file_name)
        except:
            self.debug_error('this file not could be oppened.')
        else:
            filter_seqs = type(self.target_seqs) == list and len(self.target_seqs) > 0
            filter_genes = type(self.target_genes) == list and len(self.target_genes) > 0
            lines = 0
            with open(self.file_name) as reader:
                self.genes = {}
                raw_mrnas = []
                for line in reader:
                    lines += 1
                    line = line.strip()
                    if line.startswith('#') or line.count('\t') != 8:
                        continue

                    ## To Performace: First step parse only genes
                    if not '\tgene\t' in line:
                        if '\tmRNA\t' in line:
                            raw_mrnas.append(line)
                        continue
                    
                    feature = Feature(line, style=self.style, idGen=self.idGenerator, saveRAM=saveRAM)

                    if feature.is_gene():
                        gene = Gene(line, style=self.style, idGen=self.idGenerator)

                        if filter_seqs:
                            if not gene.seqid in self.target_seqs:
                                if filter_genes and gene.id in self.target_genes:
                                    ## permite passar gene mesmo que seq nao esteja no filtro
                                    pass
                                else:
                                    continue
                        elif filter_genes and not gene.id in self.target_genes:
                            continue

                        self.genes[gene.id] = gene
                self.seqs = list(set([g.seqid for g in self.genes.values()]))
                self.debug_info('%d genes loaded in %d sequences.' % (len(self.genes), len(self.seqs)))
                
                self.mrnas = {}
                mrnas_invalids = []
                for line in raw_mrnas:
                    feature = Feature(line, style=self.style, idGen=self.idGenerator, saveRAM=saveRAM)

                    if len(feature.parents) == 1:
                        parent = feature.parents[0]
                        if parent in self.genes:
                            gene = self.genes[feature.parents[0]]
                            mrna = Mrna(line, gene, style=self.style, idGen=self.idGenerator)
                            mrna.parsed_parents[parent] = gene
                            self.mrnas[mrna.id] = mrna
                        # else: ## gene not in filter or gene is pseudo as gene-LOC113725008
                        #     self.debug_error('gene not found: ' + parent)
                    else:
                        mrnas_invalids.append(feature.id)

                if len(mrnas_invalids) > 0:
                    self.debug_warning('%d mRNAs with invalid Parent not loaded: %s' % (len(mrnas_invalids), ', '.join(mrnas_invalids)))
                    del mrnas_invalids

                self.debug_info('%s mRNAs loaded' % len(self.mrnas))
                del raw_mrnas

                self.debug_info('parsing features')
                self.features = {}
                reader.seek(0)
                cont = 0
                stat = -1
                invalids = []
                unknowns  = []
                skiped = 0
                for line in reader:
                    cont += 1
                    if lines > 10000  and int(cont / lines * 100) > stat:
                        stat = int(cont / lines * 100) 
                        self.debug_info('%d%%' % stat)
                    line = line.strip()
                    if line.startswith('#') or line.count('\t') != 8:
                        continue

                    ## To RAM Performace: Second step parse others

                    if '\tgene\t' in line or '\tmRNA\t' in line:
                        continue

                    ps = line.split('Parent=')
                    
                    if len(ps) < 2:
                        continue

                    valid_parents = [p for p in ps[-1].split(';')[0].split(',') if p in self.mrnas]

                    if len(valid_parents) < 1:
                        skiped += 1
                        continue

                    feature = Feature(line, style=self.style, idGen=self.idGenerator, saveRAM=saveRAM)
                    persist = False

                    if feature.feature in unknowns:
                        invalids.append(feature.id)
                        continue

                    for parent in valid_parents:
                        mrna = self.mrnas[parent]
                        feature.parsed_parents[parent] = mrna
                        if feature.feature == 'exon':
                            mrna.exons.append(feature)
                            persist = True
                        elif feature.feature == 'intron':
                            mrna.introns.append(feature)
                            persist = True
                        elif feature.feature == 'CDS':
                            mrna.cds.append(feature)
                            persist = True
                        elif feature.feature == 'five_prime_UTR':
                            mrna.five_prime_UTRs.append(feature)
                            persist = True
                        elif feature.feature == 'three_prime_UTR':
                            mrna.three_prime_UTRs.append(feature)
                            persist = True
                        else:
                            unknowns.append(feature.feature)
                            self.debug_warning('Unknown feature: ' + feature.feature)
                            invalids.append(feature.id)

                    if persist:
                        self.features[feature.id] = feature
                
                if skiped > 0:
                    self.debug_info('%d skiped features' % skiped)

                if len(invalids) > 0:
                    self.debug_warning('%d features invalid not loaded.' % len(invalids))
                    del invalids

                self.debug_info("%s features loaded (%d exons, %d introns, %d cds, %d 5'UTR, %d 3'UTR)" % (
                    len(self.features), 
                    len([x for x in self.features.values() if x.feature == 'exon']),
                    len([x for x in self.features.values() if x.feature == 'intron']),
                    len([x for x in self.features.values() if x.feature == 'CDS']),
                    len([x for x in self.features.values() if x.feature == 'five_prime_UTR']),
                    len([x for x in self.features.values() if x.feature == 'three_prime_UTR'])))

                    

class RefSeqGFF(GFF):
    def __init__(self, file, target_genes=None, target_seqs=None, debug=['E'], saveRAM=True):
        GFF.__init__(self, file=file, target_genes=target_genes, target_seqs=target_seqs, debug=debug, style=Style_Refseq, saveRAM=saveRAM)

