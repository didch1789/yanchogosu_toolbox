function to(varargin)

switch varargin{1}
    case 'projects'
        if ismac
            cd('/Users/jungwookim/Dropbox');
        elseif ~ismac
            cd('/media/das/cocoanlab_Dropbox/projects');
        end
    case 'dropbox'
        if ~ismac
            cd('/home/jungwoo/Dropbox')
        elseif ismac
            cd('/Users/jungwookim/Dropbox')
        end
    case 'git'
        if ~ismac 
            cd('/home/jungwoo/Dropbox/github/')
        elseif ismac
            cd('/Users/jungwookim/Dropbox/github') 
        end
    case 'data'
        if ~ismac
            cd('/media/das/cocoanlab_Dropbox/data')
        end
    otherwise
        error('No match with your input!')
end

end