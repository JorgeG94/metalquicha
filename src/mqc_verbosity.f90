!! Whether the user asked to see something once per fragment or displacement
module mqc_verbosity
   !! **One derivation, because there used to be three and they disagreed.**
   !!
   !! Output here is gated twice: by the logger level, and by a `verbose`
   !! logical threaded down call chains. The two answer different questions,
   !! and conflating them is what broke:
   !!
   !! * a **level** answers *how much* -- and it is applied where the message
   !!   is emitted, by calling `logger%large_info` rather than `logger%info`.
   !!   No predicate is needed for that, which is why there is only one
   !!   function here.
   !! * a **flag** answers *whose* calculation this is. An atomic guess, a
   !!   basis-ladder rung, a Fukui ion and a finite-difference displacement
   !!   must stay silent even at `debug`, because their output belongs to a
   !!   calculation the caller owns rather than to the one the user asked for.
   !!   A level cannot say that; those keep their flag.
   !!
   !! The one case that is genuinely both is a calculation the driver owns and
   !! runs *many* of: one SCF per fragment, one per displacement. Whether to
   !! hear those is a level question -- but it has to be answered as a flag,
   !! because it is passed into the callee. That is what this asks.
   !!
   !! It used to be asked in three places with three answers:
   !! `mqc_config_adapter` pinned the flag `.false.`, so an unfragmented run
   !! printed no method output at any level -- not its iteration table, and not
   !! `E(CASSCF)` -- while the two schedulers derived it at `verbose_level`, so
   !! per-fragment output could only be had at the loudest setting there is.
   use pic_logger, only: logger => global_logger, verbose_level
   implicit none
   private

   public :: per_item_requested

contains

   function per_item_requested() result(wanted)
      !! Whether to let a calculation the driver runs once per item speak
      !!
      !! `verbose` and above, a notch past `large_info`, because this output
      !! grows with the *count* of things rather than the size of one: a
      !! fragment SCF's iteration table on a water 20-mer is two hundred
      !! tables, and a distributed Hessian runs 6N+1 of them.
      logical :: wanted
      integer :: level

      call logger%configuration(level=level)
      wanted = level >= verbose_level
   end function per_item_requested

end module mqc_verbosity
