

/^@/{
  print;
  next;
}

$7=="="&&$1==id{
  print prev;
  print $0;
}

$7=="="{
  id=$1;
  prev=$0;
}
