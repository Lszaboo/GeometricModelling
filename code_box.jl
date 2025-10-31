
function CodeBox(name,failSafeVal)
    
    symbolName = Symbol(name)
    T = typeof(failSafeVal)
    

    txt1 = TextBox("$(name) = $(failSafeVal)")
    gd1 = GenericDependent(failSafeVal,[txt1]) do txt1
        value = nothing

        try
            text = "begin\n" * txt1[:text] * "\nend"
            textSymbols = Meta.parse(text)
        
            UserFunc = @eval function ()
                try 
                    $textSymbols
                    return $symbolName
                catch err
                    println("Error occured evaluating UserFunc:\n$(err)")
                end
                return nothing
            end

            value = Base.invokelatest(UserFunc)
        catch err
            println("Error occured parsing UserFunc:\n$(err)\nFor this text:\n$(text)")
        end

        if !(isa(value,T))
            println("$(name):\n$(value)\nWasn't a $(T)\nReturning failSafeVal...")
            return failSafeVal
        end

        return value
    end

    return gd1
end

