#!/usr/bin/python3

def blacklist():
    ''' Arguments:
        Returns:
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
{G}
{A2}
{T2}
{C2}
{G2}
'''

    As = '>ConsecutiveA\n'+'A'*50
    Ts = '>ConsecutiveT\n'+'T'*50
    Cs = '>ConsecutiveC\n'+'C'*50
    Gs = '>ConsecutiveG\n'+'G'*50
    As2 = '>ConsecutiveAb\n'+'A'*30
    Ts2 = '>ConsecutiveTb\n'+'T'*30
    Cs2 = '>ConsecutiveCb\n'+'C'*30
    Gs2 = '>ConsecutiveGb\n'+'G'*30
    Context = {
            'A': As,
            'T': Ts,
            'C': Cs,
            'G': Gs,
            'A2': As2,
            'T2': Ts2,
            'C2': Cs2,
            'G2': Gs2
            }
    goodBlacklist = Page.format(**Context)
    with open('TruSeq3-PE.fa','w') as F:
        F.write(goodBlacklist)


if __name__ == '__main__':
    blacklist()
