function to(varargin)

switch varargin{1}
    case 'projects'
        if ismac
            cd('/Users/jungwookim/Dropbox');
        elseif IsLinux
            cd('/media/das/cocoanlab_Dropbox/projects');
        end
    case 'Dropbox'
        if IsLinux
            cd('/home/jungwoo/Dropbox')
        elseif ismac
            cd('/home/jungwookim/Dropbox')
        end
    case 'git'
        if IsLinux
            cd('/home/jungwoo/github/')
        elseif ismac
            cd('/Users/jungwookim/Dropbox/github')
        end
end

end