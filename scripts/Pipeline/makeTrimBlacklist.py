#!/usr/bin/python3

def blacklist():
    ''' Arguments:
        Returns:
    '''
    Page = '''>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
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
            'G': Gs
            }
    goodBlacklist = Page.format(**Context)
    with open('TruSeq3-PE.fa','w') as F:
        F.write(goodBlacklist)


if __name__ == '__main__':
    blacklist()
