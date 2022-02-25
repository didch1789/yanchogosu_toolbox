function sendmail_ycgosu(to, subject, mssge, attach)
% work only when turning on the access of 


props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.port', '587');
props.setProperty('mail.smtp.starttls.enable', 'true');
props.setProperty('mail.smtp.auth', 'true');


setpref('Internet', 'E_mail', 'didch1789@gmail.com');
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username','didch1789');
setpref('Internet','SMTP_Password', '%%');




end