hljs.registerLanguage('shell', function(hljs) {
    var KEYWORDS = {
      keyword: 'if then else fi for while in do done exit return break continue case esac',
      literal: 'true false',
      built_in: 'echo read cd pwd mkdir hisat2 bwa'
    };
    return {
      keywords: KEYWORDS,
      contains: [
        hljs.HASH_COMMENT_MODE,
        hljs.BACKSLASH_ESCAPE
      ]
    };
  });