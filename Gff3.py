
import csv

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

    def __init__(self, file, target_genes=None, target_seqs=None, debug=['I', 'W', 'E']):
        self.file_name = file
        self.debug = debug
        
        self.target_genes = [target_genes] if type(target_genes) == str else target_genes
        self.target_seqs = [target_seqs] if type(target_seqs) == str else target_seqs

        if not self.target_genes is None:
            self.debug_info('Filtering for genes: ' + ', '.join(self.target_genes))
        
        if not self.target_seqs is None:
            self.debug_info('Filtering for seqs: ' + ', '.join(self.target_seqs))

        
        self.debug_info('Opening file {}'.format(self.file_name))

        raw_data = list(csv.reader(open(file), delimiter='\t'))
        




