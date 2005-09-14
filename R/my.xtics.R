"my.xtics" <-
function(at,length=0.5) {
  tick.y0 <- unit(0,"npc")
  tick.y1 <- unit(-length, "lines")
  grid.segments(unit(at, "native"), tick.y0,
                unit(at, "native"), tick.y1)
}

