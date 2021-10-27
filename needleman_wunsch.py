import numpy as np

def parse_blosumfile(blosum_path):
    blosum_matrix = dict()
    with open(blosum_path) as r:
        amino_acids = next(r).strip().split()
        for line in r:
            line = line.strip().split()
            amino_acid,scores = line[0],line[1:]
            for x in range(0,len(scores)):
                subscore = scores[x]
                blosum_matrix[amino_acids[x]+"_"+amino_acid] = int(subscore)
    r.close()
    return blosum_matrix

blosum_matrix = parse_blosumfile("blosum62.txt")

class global_alignment:
    def __init__(self,submatrix,gap,extendrecent,extend):
        self.submatrix = submatrix
        self.gap = gap
        self.extend = extend
        self.extendrecent = extendrecent
    def init_matrix(self):
        self.scorem  = np.zeros((self.len1,self.len2),float)
        self.scorex = np.zeros((self.len1,self.len2),float)
        self.scorey = np.zeros((self.len1,self.len2),float)
        for i in range(1,self.len1):
            self.scorem[i][0],self.scorex[i][0],self.scorey[i][0] = \
                np.NINF,np.NINF,-(self.gap+self.extend*(i-1))
        for j in range(1,self.len2):
            self.scorem[0][j],self.scorex[0][j],self.scorey[0][j] = \
                np.NINF,-(self.gap+self.extend*(j-1)),np.NINF
        self.scorex[0][0] = np.NINF
        self.scorey[0][0] = np.NINF

    def score_matrix(self,i,j,amino_key):
        X = ((-self.gap+self.scorem[i][j-1]),\
            (-self.extend+self.scorex[i][j-1]),\
            (-self.extendrecent+self.scorey[i][j-1]),\
            )
        Y = ((-self.gap+self.scorem[i-1][j]),\
            (-self.extendrecent+self.scorex[i-1][j]),\
            (-self.extend+self.scorey[i-1][j]),
        )
        M = ((self.submatrix[amino_key]+self.scorem[i-1][j-1]),\
            (self.submatrix[amino_key]+self.scorex[i-1][j-1]),\
            (self.submatrix[amino_key]+self.scorey[i-1][j-1]),\
        )
        return M,X,Y
    def score_alignment(self,peptide1,peptide2):
        self.pep1,self.pep2 = list(peptide1),list(peptide2)
        self.len1,self.len2 = len(self.pep1)+1,len(self.pep2)+1
        self.init_matrix()

    def traceback(self):
        i,j = self.len1-1,self.len2-1
        self.traceback_pep1 = list()
        self.traceback_pep2 = list()
        while i > 0 or j > 0:
            amino_key = self.pep2[j-1]+"_"+self.pep1[i-1]
            M,X,Y = self.score_matrix(i,j,amino_key)
            cell_concat = M+X+Y
            max_score = cell_concat.index(max(cell_concat))
            if max_score%3 == 0:
                #previous matrix was M
                i -= 1
                j -= 1
                self.traceback_pep2.append(self.pep2[j])
                self.traceback_pep1.append(self.pep1[i])
            elif max_score%3 == 1:
                #previous matrix was X
                j -= 1
                self.traceback_pep2.append(self.pep2[j])
                self.traceback_pep1.append("-")
            elif max_score%3 == 2:
                #previous matrix was Y
                i -= 1
                self.traceback_pep2.append("-")
                self.traceback_pep1.append(self.pep1[i])
        return "".join(self.traceback_pep1[::-1]),\
            "".join(self.traceback_pep2[::-1])
test = global_alignment(blosum_matrix,15,7,3)
test.score_alignment("MATKGTKRSYEQMETDGERQNATEIRASVGKMIDGIGRFYIQMCTELKLSDYEGRLIQNSLTIERMVLSAFDERRNKYLEEHPSAGKDPKKTGGPIYRRVDGKWRRELILYDKEEIRRIWRQANNGDDATAGLTHMMIWHSNLNDATYQRTRALVRTGMDPRMCSLMQGSTLPRRSGAAGAAVKGVGTMVMELIRMIKRGINDRNFWRGENGRRTRIAYERMCNILKGKFQTAAQRTMVDQVRESRNPGNAEFEDLIFLARSALILRGSVAHKSCLPACVYGSAVASGYDFEREGYSLVGIDPFRLLQNSQVYSLIRPNENPAHKSQLVWMACHSAAFEDLRVSSFIRGTKVVPRGKLSTRGVQIASNENMETMESSTLELRSRYWAIRTRSGGNTNQQRASSGQISIQPTFSVQRNLPFDRPTIMAAFTGNTEGRTSDMRTEIIRLMESARPEDVSFQGRGVFELSDEKATSPIVPSFDMSNEGSYFFGDNAEEYDN","MSNMDIDSINTGTIDKTPEELTPGTSGATRPIIKPATLAPPSNKRTRNPSPERTTTSSETDIGRKIQKKQTPTEIKKSVYKMVVKLGEFYNQMMVKAGLNDDMERNLIQNAQAVERILLAATDDKKTEYQKKRNARDVKEGKEEIDHNKTGGTFYKMVRDDKTIYFSPIKITFLKEEVKTMYKTTMGSDGFSGLNHIMIGHSQMNDVCFQRSKGLKRVGLDPSLISTFAGSTLPRRSGTTGVAIKGGGTLVDEAIRFIGRAMADRGLLRDIKAKTAYEKILLNLKNKCSAPQQKALVDQVIGSRNPGIADIEDLTLLARSMVVVRPSVASKVVLPISIYAKIPQLGFNTEEYSMVGYEAMALYNMATPVSILRMGDDAKDKSQLFFMSCFGAAYEDLRVLSALTGTEFKPRSALKCKGFHVPAKEQVEGMGAALMSIKLQFWAPMTRSGGNEVSGEGGSGQISCSPVFAVERPIALSKQAVRRMLSMNVEGRDADVKGNLLKMMNDSMAKKTSGNAFIGKKMFQISDKNKVNPIEIPIKQTIPNFFFGRDTAEDYDDLDY") #A,B
#test.score_alignment("MATKGTKRSYEQMETDGERQNATEIRASVGKMIDGIGRFYIQMCTELKLSDYEGRLIQNSLTIERMVLSAFDERRNKYLEEHPSAGKDPKKTGGPIYRRVDGKWRRELILYDKEEIRRIWRQANNGDDATAGLTHMMIWHSNLNDATYQRTRALVRTGMDPRMCSLMQGSTLPRRSGAAGAAVKGVGTMVMELIRMIKRGINDRNFWRGENGRRTRIAYERMCNILKGKFQTAAQRTMVDQVRESRNPGNAEFEDLIFLARSALILRGSVAHKSCLPACVYGSAVASGYDFEREGYSLVGIDPFRLLQNSQVYSLIRPNENPAHKSQLVWMACHSAAFEDLRVSSFIRGTKVVPRGKLSTRGVQIASNENMETMESSTLELRSRYWAIRTRSGGNTNQQRASSGQISIQPTFSVQRNLPFDRPTIMAAFTGNTEGRTSDMRTEIIRLMESARPEDVSFQGRGVFELSDEKATSPIVPSFDMSNEGSYFFGDNAEEYDN","MAGQGTKRTFEQMETDSKQNTTEIRSAVGRMVKAIGRFYIQMCAELKLDDKEAVLIQNSLTIERMVLSAFDERRNKYLEEHPTVGKDPKKTGGPIYRRKEGKWEREMVLMEKENIRAIWKMANDGEENLSGLSHIMIWHSNLNDTTYQRTRALVRTGMDPRMCSLMQGSTLPRRAGAAGAAIKGVGTLIMELIRMIKRGMNDRNFWKGEQGKRTRAAYERICNNLKNKFQTAPQKAMVDQVKEGKNPGNAEIEDLLFLARSALILRGAVAHKSSLPACVYGLGVSRGFDFEREGYSLVGRDPYMLLQNSQIFSIIRKGENAAHKSQLVWMACHAAAFEDIRVSSFIKGNKIVPRGKLETRGLQIAGSETLDEALVVSLDIKSHYWAIKTRSGGNPQQSRSSAGQIAVQPTFSVQRNIPFEKKTIMAAFSNIEEGRITDMRTEIIKLMENSDPKDKVFLGRGVFEMADEKATNPIVPSLDGNDEGSYFFGDKAEEFDI") #A,BatA
#test.score_alignment("MATKGTKRSYEQMETDGERQNATEIRASVGKMIDGIGRFYIQMCTELKLSDYEGRLIQNSLTIERMVLSAFDERRNKYLEEHPSAGKDPKKTGGPIYRRVDGKWRRELILYDKEEIRRIWRQANNGDDATAGLTHMMIWHSNLNDATYQRTRALVRTGMDPRMCSLMQGSTLPRRSGAAGAAVKGVGTMVMELIRMIKRGINDRNFWRGENGRRTRIAYERMCNILKGKFQTAAQRTMVDQVRESRNPGNAEFEDLIFLARSALILRGSVAHKSCLPACVYGSAVASGYDFEREGYSLVGIDPFRLLQNSQVYSLIRPNENPAHKSQLVWMACHSAAFEDLRVSSFIRGTKVVPRGKLSTRGVQIASNENMETMESSTLELRSRYWAIRTRSGGNTNQQRASSGQISIQPTFSVQRNLPFDRPTIMAAFTGNTEGRTSDMRTEIIRLMESARPEDVSFQGRGVFELSDEKATSPIVPSFDMSNEGSYFFGDNAEEYDN","MASQGTKRSYEQMETGGERQNATEIRASVGRMVGGIGRFYIQMCTELKLSDYEGRLIQNSITIERMVLSAFDERRNKYLEEHPSAGKDPKKTGGPIYRRRDGKWVRELILYDKEEIRRIWRQANNGEDATAGLTHMMIWHSNLNDATYQRTRALVRTGMDPRMCSLMQGSTLPRRSGAAGAAVKGVGTMVMELIRMIKRGINDRNFWRGENGRRTRIAYERMCNILKGKFQTAAQRAMMDQVRESRNPGNAEIEDLIFLARSALILRGSVAHKSCLPACVYGLAVASGYDFEREGYSLVGIDPFRLLQNSQVFSLIRPNENPAHKSQLVWMACHSAAFEDLRVSSFIRGTRVAPRGQLSTRGVQIASNENMETMDSSTLELRSRYWAIRTRSGGNTNQQRASAGQISVQPTFSVQRNLPFERATIMAAFTGNTEGRTSDMRTEIIRMMESSRPEDVSFQGRGVFELSDEKATNPIVPSFDMSNEGSYFFGDNAEEYDN") #A,H5N1A

print test.traceback()
