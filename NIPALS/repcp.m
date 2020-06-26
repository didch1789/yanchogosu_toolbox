function repcp(out, in)
    % out: string you want to change. you can also enter 'digit' for
    % numbers, 'wspace' for blanks, 'lower' for lower case letter, 'upper'
    % for upper case letter
    %
    % in: string you want to replace
    
    copy_text = char(clipboard('paste'));
    if contains(out, {'digit', 'wspace', 'lower', 'upper'})
        idx = isstrprop(copy_text, out);
        copy_text(idx) = in;
        paste_text = copy_text;
    else
        paste_text = strrep(copy_text, out, in);     
    end
    clipboard('copy', paste_text)
end