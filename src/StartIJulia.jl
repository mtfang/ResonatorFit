# Short script which starts an IJulia notebook in the folder "<User>\IJulia"
# The folder "<User>\IJulia" will be created if it does not exist

function installedIJulia()
    if typeof(Pkg.installed("IJulia")) == VersionNumber return true
    else return false
    end
end

function updateRebuild()
    Pkg.update()
    Pkg.build("IJulia")
end

if !installedIJulia()
    Pkg.add("IJulia")
end

using IJulia
# On Windows $(homedir()) corresponds to C:\Users\<User>
mkpath("$(homedir())/IJulia")
cd("$(homedir())/IJulia")
notebook()
