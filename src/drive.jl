using ProgressMeter

function isfile_drive(parentId, fileName)
    length(readall(`drive list -n -q "'$parentId' in parents and title = '$fileName'"`)) != 0
end

function dir2drive(sourceDir, destName)
    destId = split(readall(`drive list -n --title $destName`))[1]
    @showprogress for fileName in readdir(sourceDir)
        if !isfile_drive(destId, fileName)
            readall(`drive upload -f $sourceDir/$fileName --title $fileName --parent $destId`)
        end
    end
end
