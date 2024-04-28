function st = randomseq_generator(strlen)

   symbols = ['a':'z' 'A':'Z' '0':'9'];
   nums = randi(numel(symbols),[1 strlen]);
   st = symbols (nums);

end
