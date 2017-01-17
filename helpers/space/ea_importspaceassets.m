function ea_importspaceassets(~,~,fromspace,what)

ea_genwarp2space(fromspace);


switch what
    case 'both'
        ea_warplabelassets(fromspace)
        ea_warpatlassets(fromspace)
    case 'atlases'
        ea_warpatlassets(fromspace)
        
    case 'labeling'
        ea_warplabelassets(fromspace)
        
        
end


function ea_warplabelassets(fromspace)

keyboard

function ea_warpatlasassets(fromspace)



