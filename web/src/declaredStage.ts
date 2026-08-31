/**
 * What the operator says the material is, set against what the array says it is.
 *
 * THE OPERATOR KNOWS WHAT THEY PUT IN THE MACHINE. The stage inference reads heterozygosity, call
 * rate and allele-frequency spread off the array and guesses from those, which is the only thing it
 * can do with no other information, and it is wrong often enough to matter: an array that failed its
 * reaction imitates haploid material, and a genuinely haploid product with drop-in imitates a
 * diploid one. Somebody who dissected the cell does not have to guess.
 *
 * SO THE DECLARATION WINS, AND THE DISAGREEMENT IS REPORTED. Overriding silently would let a
 * mislabelled tube produce a confident answer with nothing to warn the reader. Refusing to be
 * overridden would make the tool argue with the person holding the pipette. Reporting both is the
 * only version that is useful and honest at once.
 *
 * THE PLOIDY AXIS IS THE ONE THAT DECIDES ANYTHING. Which named stage a diploid cell is changes
 * detection floors. Whether the cell is diploid at all changes whether parent of origin is a
 * question that can be asked: a genome carrying one parental contribution has no second parent to
 * attribute anything to, and every channel that names a parent by comparing the two is meaningless
 * on it. A haploid-versus-diploid disagreement is therefore not a detail, and it is reported
 * separately from a disagreement between two diploid stages.
 */
import type { Stage, StageCall } from './stage.ts'

/** Whether a stage names material with one parental genome in it or two. */
export type Ploidy = 'haploid' | 'diploid' | 'unknown'

/**
 * A stage's ploidy.
 *
 * `single-cell` is an amplified single cell of a diploid line, so it is diploid however much
 * dropout it carries. `failed` and `unknown` name the absence of a call, not a genome.
 */
export function ploidyOf(stage: Stage): Ploidy {
  switch (stage) {
    case 'bulk':
    case 'trophectoderm':
    case 'blastomere':
    case 'single-cell':
      return 'diploid'
    case 'haploid':
      return 'haploid'
    case 'failed':
    case 'unknown':
      return 'unknown'
    default: {
      const unhandled: never = stage
      return unhandled
    }
  }
}

export interface StageAgreement {
  /** The stage the run will actually use. The declaration where there is one. */
  used: StageCall
  /** What the array reads, always computed, even when it is overridden. */
  inferred: StageCall
  declared?: Stage
  agrees: boolean
  /** Set when the two differ on whether the genome carries one parental contribution or two. */
  ploidyConflict: boolean
  /** Empty when they agree. Written for the operator, not for a log. */
  notice: string
}

/**
 * A declared stage, carrying the inferred stage's measurements where they still apply.
 *
 * Dropout comes from the DECLARED stage, because dropout is a property of how much template the
 * material had, which is what the declaration is about. Keeping the inferred array's dropout under
 * a declaration would parameterise every likelihood below with a figure for material the operator
 * has just said this is not.
 */
function asDeclared(declared: Stage, inferred: StageCall, defaults: StageCall): StageCall {
  return {
    ...defaults,
    stage: declared,
    basis: 'stage-default',
    why: `Declared ${declared} by the operator. The array reads ${inferred.stage}: ${inferred.why}`,
    caveat: inferred.caveat,
  }
}

/**
 * Reconcile what the operator declared with what the array reads.
 *
 * `defaultsFor` supplies the shipped StageCall for a named stage, so this module does not carry a
 * second copy of the dropout table.
 */
export function reconcileStage(
  inferred: StageCall,
  declared: Stage | undefined,
  defaultsFor: (s: Stage) => StageCall,
): StageAgreement {
  if (!declared) {
    return { used: inferred, inferred, agrees: true, ploidyConflict: false, notice: '' }
  }
  if (declared === inferred.stage) {
    return {
      used: inferred, inferred, declared, agrees: true, ploidyConflict: false,
      notice: `Declared ${declared}, and the array reads the same.`,
    }
  }

  const dp = ploidyOf(declared)
  const ip = ploidyOf(inferred.stage)
  const ploidyConflict = dp !== ip && dp !== 'unknown' && ip !== 'unknown'

  const head = `Declared ${declared}, but this array reads ${inferred.stage}. `
    + `The declaration is being used. ${inferred.why}`

  // A disagreement between two diploid stages moves detection floors. A disagreement about whether
  // the genome is diploid at all moves whether parent of origin is answerable, so it is said louder
  // and in both directions, because which of the two is wrong changes what the reader should do.
  const notice = !ploidyConflict
    ? `${head} Both name material with two parental genomes in it, so parent of origin remains `
      + 'answerable either way; what changes is the detection floor each event is measured against.'
    : dp === 'diploid'
      ? `${head} THESE DISAGREE ABOUT PLOIDY. You have said this cell carries both parents; the `
        + 'array reads one parental genome only. If the declaration is right, the array is '
        + 'under-calling heterozygosity, usually amplification dropout, and parental origin is '
        + 'answerable but every heterozygous-site statistic on it is depressed. If the array is '
        + 'right, the sample is not what the tube says and no comparison between two parental '
        + 'genomes in it means anything. Resolve this before using any parental call from this row.'
      : `${head} THESE DISAGREE ABOUT PLOIDY. You have said this carries one parental genome; the `
        + 'array reads two. If the declaration is right, the extra heterozygosity is drop-in or a '
        + 'second cell in the tube, and origin calls made by comparing two parental genomes here '
        + 'would be comparing a contaminant. If the array is right, this is biparental material '
        + 'and more is answerable on it than the declaration allows.'

  return {
    used: asDeclared(declared, inferred, defaultsFor(declared)),
    inferred,
    declared,
    agrees: false,
    ploidyConflict,
    notice,
  }
}
