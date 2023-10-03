  update <- function() {
    isolate({
      # get HEADERS/ANNOTS (raw epochs)
      ret <- leval("SEGMENTS & HEADERS & ANNOTS & DUMP-MASK")

      # no records left?
      has.records <- "SEGMENTS" %in% names(ret)
      if (has.records && !any(ret$DUMP_MASK$E$EMASK == 0)) has.records <- F

      # check records set?
      if (!has.records) {
        showModal(modalDialog(
          title = "No unmasked records left",
          "Please Refresh or reload a valid EDF",
          easyClose = TRUE
        ))
      }

      req(has.records == T)

      values$opt[["header1"]] <- ret$HEADERS$BL
      values$opt[["header1"]]$EPOCH <- values$opt[["header1"]]$REC_DUR_SEC / values$elen

      values$opt[["header2"]] <- ret$HEADERS$CH
      values$opt[["header2"]] <- values$opt[["header2"]][, c("CH", "PDIM", "SR", "PMIN", "PMAX", "TRANS")]
      values$opt[["header2"]]$PMIN <- signif(values$opt[["header2"]]$PMIN, 4)
      values$opt[["header2"]]$PMAX <- signif(values$opt[["header2"]]$PMAX, 4)
      names(values$opt[["header2"]]) <- c("Channel", "Unit", "SRate", "Min", "Max", "Transducer")

      # segments
      values$opt[["curr.segsumm"]] <- ret$SEGMENTS$BL
      values$opt[["curr.segidx"]] <- ret$SEGMENTS$SEG[, c("START", "STOP")]
      values$opt[["curr.segments"]] <- ret$SEGMENTS$SEG[, c("SEG", "START_HMS", "STOP_HMS", "DUR_MIN", "DUR_SEC", "START", "STOP")]
      values$opt[["curr.segments"]]$SEG <- paste("Seg", values$opt[["curr.segments"]]$SEG)

      if (ret$SEGMENTS$BL$NGAPS > 0) {
        t2 <- ret$SEGMENTS$GAP[, c("GAP", "START_HMS", "STOP_HMS", "DUR_MIN", "DUR_SEC", "START", "STOP")]
        t2$GAP <- paste("Gap", t2$GAP, sep = " ")
        names(t2)[1] <- "SEG"
        values$opt[["curr.segments"]] <- rbind(values$opt[["curr.segments"]], t2)
        values$opt[["curr.segments"]] <- values$opt[["curr.segments"]][order(values$opt[["curr.segments"]]$START), ]
      }

      values$opt[["curr.segments"]]$START <- round(values$opt[["curr.segments"]]$START, 2)
      values$opt[["curr.segments"]]$STOP <- round(values$opt[["curr.segments"]]$STOP, 2)
      values$opt[["curr.segments"]]$DUR_SEC <- round(values$opt[["curr.segments"]]$DUR_SEC, 2)
      values$opt[["curr.segments"]]$DUR_MIN <- round(values$opt[["curr.segments"]]$DUR_MIN, 2)
      names(values$opt[["curr.segments"]]) <- c("Segment", "Start", "Stop", "Duration (m)", "Duration (s)", "Start(s)", "Stop(s)")

      # Channels
      values$opt[["chs"]] <- ret$HEADERS$CH$CH

      # Sample rates
      values$opt[["sr"]] <- ret$HEADERS$CH$SR
      names(values$opt[["sr"]]) <- as.character(ret$HEADERS$CH$CH)

      # Channel units
      values$opt[["units"]] <- ret$HEADERS$CH$PDIM
      names(values$opt[["units"]]) <- as.character(ret$HEADERS$CH$CH)

      # Channel type
      values$opt[["type"]] <- ret$HEADERS$CH$TYPE

      # Annots
      #  (skip 'Sleep Stage' and special annots)
      skips <- c("SleepStage", "duration_hms", "duration_sec", "epoch_sec", "start_hms")
      values$opt[["annots"]] <- ret$ANNOTS$ANNOT$ANNOT[!ret$ANNOTS$ANNOT$ANNOT %in% skips]
      values$hasannots <- length(values$opt[["annots"]]) > 0
      values$opt[["annots.summ"]] <- ret$ANNOTS$ANNOT[ret$ANNOTS$ANNOT$ANNOT != "SleepStage", ]
      values$opt[["annots.inst"]] <- ret$ANNOTS$ANNOT_INST_T1_T2[ret$ANNOTS$ANNOT_INST_T1_T2$ANNOT != "SleepStage", ]

      # Get mask
      values$opt[["unmasked"]] <- ret$DUMP_MASK$E$E[ret$DUMP_MASK$E$EMASK == 0]
      values$opt[["included"]] <- ret$DUMP_MASK$E$E

      # clear misc other
      values$estats <- NULL
    })
  }
