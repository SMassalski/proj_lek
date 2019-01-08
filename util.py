
import requests



class Uniprot:
    BASE = 'http://www.uniprot.org'
    KB_ENDPOINT = '/uniprot/'
    TOOL_ENDPOINT = '/uploadlists/'
    #cols: separated by commas, no space
    def query(q_str,fmt='list',cols=None,random='no'):
        params = {'query':q_str,'format':fmt,'random':random}
        if fmt=='tab' and cols:
            params['columns'] = cols
        return requests.get(Uniprot.BASE + Uniprot.KB_ENDPOINT,params=params)
    
    def retrieve(UPIDs,source_fmt='ACC+ID', output_fmt='fasta',cols=None):
        ID_str = ' '.join(UPIDs)
        params = {'query':ID_str,
                 'from':source_fmt,
                 'to':'ACC',
                 'format':output_fmt}
        if output_fmt=='tab' and cols:
            params['columns'] = cols
        
        return requests.get(Uniprot.BASE+Uniprot.TOOL_ENDPOINT,params=params)
    
    def parse(data):
        pass
        

