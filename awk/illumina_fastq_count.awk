BEGIN{
  filter = 0;
}
{ 
  if ( FNR % 4 == 1 && substr($2, 3, 1) == "Y" ) {
    filter+=1
  }
}
END{
  TOTAL = NR / 4;
  if ( paired == 1 ) {
    TOTAL *= 2;
    filter *= 2;
  }
  print "total", TOTAL ;
  print "pf", TOTAL - filter ;
  print "qc", filter ;
}
