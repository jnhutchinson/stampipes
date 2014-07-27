BEGIN{
  OFS="\t";
  chr=""
}
{ 
  if( $4 != 0 && $2 != 0 ){
    if( chr=="" ){
      chr=$1;
      print "variableStep chrom=" chr " span=" binI;
    }
    if( $1==chr ){
      print $2,$4;
    } else {
      chr=$1;
      print "variableStep chrom=" chr " span=" binI;
      print $2, $4;
    }
  }
}
