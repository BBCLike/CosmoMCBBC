module snovae
    implicit none
    contains
    subroutine SNLikelihood_Add(LikeList, Ini)
        use SNLS
        use Union2
        use likelihood
        use settings
        use JLA
        use BBC
        class(TLikelihoodList) :: LikeList
        class(TSettingIni) :: ini
        integer count

        if (.not. Ini%Read_Logical('use_SN',.false.)) then
            return
        endif
        count = LikeList%Count
        call SNLSLikelihood_Add(LikeList, Ini)
        call JLALikelihood_Add(LikeList, Ini)
        call BBCLikelihood_Add(LikeList, Ini)
        call Union2Likelihood_Add(LikeList, Ini)
        if (LikeList%Count > count+1) then
            call MpiStop('SNLikelihood_Add: more than one - datasets not independent')
        endif
    end subroutine SNLikelihood_Add
end module snovae
