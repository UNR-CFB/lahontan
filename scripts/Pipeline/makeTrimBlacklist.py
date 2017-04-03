#!/usr/bin/python3

'''Usage:
    makeTrimBlacklist.py [-t <newadapterfile>]

Used to create Trimmomatic blacklist called TruSeq3-PE.fa by default.
Name and path can be specified by -t or --tofile

Options:
    -h, --help
        Show this screen and exit
    -t <newadapterfile>
        Where to save new blacklist [default: ./TruSeq3-PE.fa]
'''

from docopt import docopt

def blacklist(filePath):
    ''' Arguments:
            None
        Returns:
            None

        Writes new blacklist to ./TruSeq3-PE.fa
        To edit blacklist, simply edit Page variable
        making sure to keep consistent syntax
    '''
    Page = '''>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PE1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>PE2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE2_rc
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>ABISolid3AdapterB
CCTATCCCCTGTGTGCCTTGGCAGTCTCAGCCTCTCTATGGGCAGTCGGT
>PCR_Primer1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PCR_Primer1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
>PCR_Primer2
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>PCR_Primer2_rc
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
>FlowCell1
TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC
>FlowCell2
TTTTTTTTTTCAAGCAGAAGACGGCATACGA
{A}
{T}
{C}
{G}'''
    As = '>ConsecutiveA\n'+'A'*500
    Ts = '>ConsecutiveT\n'+'T'*500
    Cs = '>ConsecutiveC\n'+'C'*500
    Gs = '>ConsecutiveG\n'+'G'*500
    Context = {
            'A': As,
            'T': Ts,
            'C': Cs,
            'G': Gs,
            }
    goodBlacklist = Page.format(**Context)
    with open(filePath,'w') as F:
        F.write(goodBlacklist)


if __name__ == '__main__':
    argument = docopt(__doc__, version='makeTrimBlacklist.py 1.00')
    blacklist(argument['-t'])
