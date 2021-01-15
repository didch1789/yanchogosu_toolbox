function to(varargin)

switch varargin{1}
    case 'projects'
        if ismac
            cd('/Users/jungwookim/Dropbox');
        elseif IsLinux
            cd('/media/das/cocoanlab_Dropbox/projects');
        end
    case 'dropbox'
        if IsLinux
            cd('/home/jungwoo/Dropbox')
        elseif ismac
            cd('/Users/jungwookim/Dropbox')
        end
    case 'git'
        if IsLinux
            cd('/home/jungwoo/github/')
        elseif ismac
            cd('/Users/jungwookim/Dropbox/github') 
        end
    case 'data'
        if IsLinux
            cd('/media/das/cocoanlab_Dropbox/data')
        end
    otherwise
        error('No match with your input!')
end

end